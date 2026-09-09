"""AreTomo3 run -> CETS ``Region`` (+ companion entries) for one series."""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from cryoet_alignment.io.cets import ctf as cets_ctf
from cryoet_alignment.io.cets.alignment import ReferenceVolume, alignment_to_cets
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import (
    AlignmentCompanion,
    ImageCompanion,
    TiltSeriesCompanion,
    TomogramCompanion,
)
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.entities import (
    movie_stack_series_entity,
    region_entity,
    tilt_series_entity,
    tomogram_entity,
)
from cryoet_alignment.io.cets.frames import FRAME_CONVENTIONS, image_frame
from cryoet_alignment.io.cryoet_data_portal import Alignment

from cets_aretomo3.discover import AreTomo3Run
from cets_aretomo3.dose import constant_exposure, pre_exposure

ALIGNMENT_NAME = "aretomo3"


@dataclass
class SeriesResult:
    region: Any
    tilt_series_companion: TiltSeriesCompanion
    alignment_companion: AlignmentCompanion
    tomogram_companions: Dict[str, TomogramCompanion]


def _rel(path: Optional[Path], base: Path, mode: str) -> Optional[str]:
    if path is None:
        return None
    p = Path(path)
    if mode == "absolute":
        return str(p.resolve())
    return os.path.relpath(p.resolve(), base.resolve())


def aretomo3_to_cets(run: AreTomo3Run, res: Resolver, sr: SeriesReport, *, out_dir: Path) -> SeriesResult:
    """Build the CETS region for one discovered AreTomo3 series. ``res`` resolves every value (CLI > config >
    discovered > default-with-warning); ``sr`` collects provenance, gates and dropped items."""
    aln = run.aln
    stem = run.stem
    n_raw = aln.n_raw

    pix = res.require("pix", discovered=run.get("pix"), note=run.source("pix") or "")
    paths_mode = res.value("paths", default="relative")
    drop_locals = bool(res.optional("drop_locals", absent=False))
    no_ctf = bool(res.optional("no_ctf", absent=False))
    dose_convention = res.value("dose_convention", default="exclusive")

    if not aln.is_rigid:
        msg = f"{aln.NumPatches} patches x {len(aln.GlobalAlignments)} tilts of local alignment"
        if not drop_locals:
            raise ValueError(f"{stem}: the .aln carries a local alignment ({msg}); the rigid profile cannot represent it - pass --drop-locals to keep the global rows only")
        sr.dropped.append(f"local alignment ({msg})")

    # --- tilt images ---------------------------------------------------------------------------
    width, height = run.get("image_dims_px")
    stage = [float(v) for v in res.value("stage_tilt_deg", discovered=run.get("stage_tilt_deg"), note=run.source("stage_tilt_deg") or "")]
    acq = res.optional("acq_index_1b", discovered=run.get("acq_index_1b"), note=run.source("acq_index_1b") or "")
    exposure = run.get("exposure")
    dose_per_tilt = res.optional("dose_per_tilt")
    if dose_per_tilt is not None:
        exposure = constant_exposure(n_raw, float(dose_per_tilt))
    elif exposure is not None:
        res.resolve("exposure", discovered=exposure, note=run.source("exposure") or "")
    pre = None
    if acq is not None and exposure is not None:
        pre = pre_exposure([int(a) for a in acq], [float(e) for e in exposure], dose_convention)
    else:
        sr.warnings.append("no acquisition order + per-image exposure (from _TLT.txt / mdoc / --dose-per-tilt): accumulated_dose left null")

    ctfs: Optional[List[Any]] = None
    df_hand = None
    if not no_ctf and run.get("ctf") is not None:
        c = run.get("ctf")
        res.resolve("ctf", discovered=run.ctf_path.name, note="_CTF.txt")
        ctfs = [cets_ctf.from_aretomo3_row(r) for r in c.rows]
        if c.has_df_hand:
            hands = {r.df_hand for r in c.rows}
            df_hand = hands.pop() if len(hands) == 1 else None

    frame_names = run.get("frame_names")
    movie_ids = [f"{stem}_movie_{i}" for i in range(n_raw)] if frame_names else None
    stack_rel = _rel(run.stack_path, out_dir, paths_mode) if run.stack_path else None
    ts = tilt_series_entity(
        tilt_series_id=stem,
        path=stack_rel,
        width=width,
        height=height,
        pixel_size_a=pix,
        nominal_angles=stage,
        doses=pre,
        ctfs=ctfs,
        movie_stack_ids=movie_ids,
        movie_stack_series_id=f"{stem}_movies" if frame_names else None,
    )
    movie_series = []
    if frame_names:
        movie_series.append(movie_stack_series_entity(
            series_id=f"{stem}_movies",
            stacks=[{"id": movie_ids[i], "path": frame_names[i]} for i in range(n_raw)],
        ))

    # --- reference volume: the alignment's own bin-1 box (always) + the reconstructed file when declared ---
    vol_z = res.require(
        "tomo_size", discovered=(run.get("tomo_dims_px")[2] * run.get("bin") if run.get("tomo_dims_px") else None),
        note=(run.source("tomo_dims_px") + " x bin" if run.source("tomo_dims_px") else ""),
    )
    vol_z_px = int(round(float(vol_z)))
    hub = Alignment.from_aretomo3(aln, vol_size_px=(width, height, vol_z_px), pixel_size_a=pix)
    ref_tomo = tomogram_entity(
        tomogram_id=f"{stem}_volume", path=None, size_px=(width, height, vol_z_px), voxel_size_a=pix, tilt_series_id=stem,
    )
    tomograms = [ref_tomo]
    tomo_companions = {
        ref_tomo.id: TomogramCompanion(voxel_implied_a=pix, source_ref="alignment box (RawSize x pixel, Z from --tomo-size or _Vol.mrc x bin); no file"),
    }
    if run.vol_path is not None and run.get("vol_layout") == "xyz":
        flip_vol = res.optional("flip_vol")
        if flip_vol in (1, 2, "1", "2"):
            dims = run.get("tomo_dims_px")
            header_voxel = run.get("vol_voxel_header_a")
            implied = pix * width / dims[0]
            if not header_voxel or header_voxel <= 0:
                sr.warnings.append(f"{run.vol_path.name}: no voxel size in the header; using the raw-extent value {implied:.5f} Å")
                header_voxel = implied
            elif abs(header_voxel - implied) > 1e-3 * implied:
                sr.warnings.append(f"{run.vol_path.name}: header voxel {header_voxel:.5f} Å differs from the raw-extent value {implied:.5f} Å; the header value is used")
            # the declared (header) voxel size is used; the raw-extent value goes to the companion for information
            tomo = tomogram_entity(
                tomogram_id=f"{stem}_tomo", path=_rel(run.vol_path, out_dir, paths_mode), size_px=dims, voxel_size_a=float(header_voxel),
                tilt_series_id=stem,
            )
            tomograms.append(tomo)
            tomo_companions[tomo.id] = TomogramCompanion(
                voxel_header_a=run.get("vol_voxel_header_a"), voxel_implied_a=implied, flip_vol=int(flip_vol),
                reconstruction_software="AreTomo3", source_ref=str(run.vol_path.name),
            )
            if int(flip_vol) == 2:
                sr.warnings.append("-FlipVol 2 volumes are axis-reversed relative to -FlipVol 1; the grid offset/hand is recorded in the companion, not modelled")
        else:
            sr.warnings.append(f"{run.vol_path.name} found but -FlipVol not declared (--flip-vol 1|2): no Tomogram entity for the file")
    elif run.vol_path is not None:
        sr.warnings.append(f"{run.vol_path.name} is an XZY (-FlipVol 0) volume: not described (declare --flip-vol and use an XYZ volume)")

    reference = ReferenceVolume.from_tomogram(ref_tomo)
    cets_alignment = alignment_to_cets(
        hub, tilt_series_id=stem, alignment_name=ALIGNMENT_NAME, image=image_frame(ts.images[0]), reference=reference,
        frame=FRAME_CONVENTIONS["ARETOMO3"],
    )
    sr.gates.append(Gate("rows", len(cets_alignment.projection_alignments) == len(aln.GlobalAlignments),
                         value=len(cets_alignment.projection_alignments), expected=len(aln.GlobalAlignments)))
    sr.gates.append(Gate("dark_sections", True, value=aln.dark_indices()))

    # --- companion --------------------------------------------------------------------------------
    images = {}
    for i in range(n_raw):
        images[ts.images[i].id] = ImageCompanion(
            acquisition_index_1b=None if acq is None else int(acq[i]),
            exposure_dose=None if exposure is None else float(exposure[i]),
            stage_angle_deg=stage[i],
            frame_name=None if frame_names is None else frame_names[i],
            use_tilt=i not in set(aln.dark_indices()),
        )
    ts_comp = TiltSeriesCompanion(
        source_tool="AreTomo3",
        voltage_kv=res.optional("voltage", discovered=run.get("voltage"), note=run.source("voltage") or ""),
        cs_mm=res.optional("cs"),
        amplitude_contrast=res.optional("amp_contrast"),
        pixel_size_acquisition_a=pix,
        tilt_axis_nominal_deg=run.get("tilt_axis_nominal_deg"),
        alpha_offset_deg=run.get("alpha_offset_deg"),
        beta_offset_deg=run.get("beta_offset_deg"),
        defocus_hand=df_hand,
        defocus_hand_convention="aretomo3 _CTF.txt dfHand" if df_hand is not None else None,
        collection_metadata_path=_rel(run.mdoc_path, out_dir, paths_mode) if run.mdoc_path else None,
        images=images,
    )
    dropped = list(sr.dropped)
    if run.get("beta_offset_deg"):
        dropped.append(f"BetaOffset {run.get('beta_offset_deg')} (CTF-only in AreTomo3; recorded)")
    aln_comp = AlignmentCompanion(
        name=ALIGNMENT_NAME, tilt_series_id=stem, format="ARETOMO3",
        alignment_type="GLOBAL", is_portal_standard=True,
        reference_tomogram_id=ref_tomo.id, tomogram_ids=[t.id for t in tomograms],
        native_volume_dimension_a=hub.volume_dimension,
        frame_convention={"image_center": "half", "volume_center": "half"}, dropped=dropped,
        thickness_px=run.get("thickness_px"), source_ref=str(run.aln_path.name),
    )
    region = region_entity(region_id=stem, tilt_series=[ts], alignments=[cets_alignment], tomograms=tomograms,
                           movie_stack_series=movie_series)
    sr.provenance = res.provenance()
    return SeriesResult(region, ts_comp, aln_comp, tomo_companions)
