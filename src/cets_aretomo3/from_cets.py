"""CETS region -> AreTomo3 ``-Cmd 2`` input directory (``<stem>.aln``, ``<stem>_TLT.txt``, ``<stem>_CTF.txt``)."""

import os
from pathlib import Path
from typing import Dict, Optional

from cryoet_alignment.io.aretomo3 import AreTomo3CTF, AreTomo3TLT
from cryoet_alignment.io.cets import ctf as cets_ctf
from cryoet_alignment.io.cets.alignment import (
    ReferenceVolume,
    alignment_from_cets,
    alignment_name_of,
    select_alignment,
    select_tomogram,
)
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.frames import FRAME_CONVENTIONS, find_by_id, image_frame

ALN_HEADER = "# AreTomo Alignment / Priims bprmMn"


def aretomo3_hint(
    root: Path,
    stem: str,
    *,
    pixel_size_a=None,
    voltage_kv=None,
    cs_mm=None,
    amplitude_contrast=None,
    vol_z_px=None,
    has_ctf: bool = False,
) -> str:
    """The ``-Cmd 2`` line for the written directory (arewarpion ``project/aretomo.py::aretomo3_hint``).
    ``-AtBin`` is the user's choice; ``-CorrCTF 1`` is suggested only when a ``_CTF.txt`` was written."""

    def f(v, spec):
        return f"{v:{spec}}" if v is not None else "<?>"

    return (
        f"AreTomo3 -Cmd 2 -Serial 0 -InPrefix {root / stem}.mrc -OutDir {root / 'out'}/ -Gpu 0 "
        f"-PixSize {f(pixel_size_a, 'g')} -kV {f(voltage_kv, 'g')} -Cs {f(cs_mm, 'g')} "
        f"-AmpContrast {f(amplitude_contrast, 'g')} -VolZ {f(vol_z_px, 'd') if vol_z_px is not None else '<?>'} "
        f"-AtBin <bin> -FlipVol 1 -Wbp 1 -CorrCTF {1 if has_ctf else 0} -SplitSum 0"
    )


def _resolve_doc_path(p: Optional[str], doc_dir: Path) -> Optional[Path]:
    if not p:
        return None
    q = Path(p)
    return q if q.is_absolute() else (doc_dir / q)


def cets_to_aretomo3(
    region,
    res: Resolver,
    sr: SeriesReport,
    *,
    out_dir: Path,
    doc_dir: Path,
    companion: Optional[Companion],
    alignment_selector=None,
    tomogram_selector: Optional[str] = None,
    overwrite: bool = False,
) -> Dict[str, Path]:
    cets_alignment = select_alignment(region, alignment_selector)
    ts = find_by_id(region.tilt_series, cets_alignment.tilt_series_id, "tilt series")
    stem = ts.id
    aln_name = alignment_name_of(cets_alignment)
    comp_aln = companion.alignment(ts.id, aln_name) if companion else None
    comp_ts = companion.tilt_series.get(ts.id) if companion else None

    tomo = select_tomogram(region, cets_alignment, tomogram_selector, companion)
    reference = ReferenceVolume.from_tomogram(tomo)
    native = comp_aln.native_volume_dimension_a if comp_aln and comp_aln.native_volume_dimension_a else None
    hub = alignment_from_cets(
        cets_alignment,
        tilt_series=ts,
        reference=reference,
        target_frame=FRAME_CONVENTIONS["ARETOMO3"],
        native_dimension_a=native,
        format_="ARETOMO3",
    )
    res.resolve(
        "reference_tomogram", discovered=tomo.id, note="companion" if comp_aln and comp_aln.tomogram_ids else "region"
    )

    # representability: X rotation
    if hub.has_x_rotation:
        worst = max(abs(p.volume_x_rotation) for p in hub.per_section_alignment_parameters)
        if not res.optional("drop_x_rotation", absent=False):
            raise ValueError(
                f"{stem}: the alignment carries a volume X rotation (max |xrot| = {worst:.4f} deg) that a .aln cannot "
                "represent; pass --drop-x-rotation to zero it explicitly",
            )
        for p in hub.per_section_alignment_parameters:
            p.volume_x_rotation = 0.0
        sr.dropped.append(f"volume X rotation up to {worst:.4f} deg")
        sr.gates.append(Gate("x_rotation", True, value=worst, note="dropped on request"))
    else:
        sr.gates.append(Gate("x_rotation", True, value=0.0))

    images = sorted(ts.images or [], key=lambda im: im.section)
    n_raw = len(images)
    if [im.section for im in images] != list(range(n_raw)):
        raise ValueError(f"{stem}: tilt image sections are not 0..{n_raw - 1}")
    fr = image_frame(images[0])
    pix = fr.isotropic_spacing
    width, height = fr.size_px
    aligned = {p.z_index for p in hub.per_section_alignment_parameters}
    dark_angles = {im.section: float(im.nominal_tilt_angle or 0.0) for im in images if im.section not in aligned}
    for im in images:
        if im.section not in aligned and im.nominal_tilt_angle is None:
            sr.warnings.append(f"dark image {im.id}: no nominal_tilt_angle, DarkFrame angle written as 0")

    z_px = int(round(reference.extent_a[2] / pix))
    vol_z = res.value("tomo_size", discovered=z_px, note=f"{tomo.id} extent / pixel")
    thickness = comp_aln.thickness_px if comp_aln else None  # AreTomo3's estimated sample thickness, never the box
    aln = hub.to_aretomo(ts_size=(width, height, n_raw), dark_angles=dark_angles, thickness_px=thickness)
    aln.header = ALN_HEADER

    # _TLT.txt: acquisition order + exposure from companion / CLI
    acq = None
    exposure = None
    if comp_ts and comp_ts.images:
        acqs = [comp_ts.images.get(im.id) for im in images]
        if all(a is not None and a.acquisition_index_1b is not None for a in acqs):
            acq = [a.acquisition_index_1b for a in acqs]
            res.resolve("acq_order", companion=acq, note="companion")
        if all(a is not None and a.exposure_dose is not None for a in acqs):
            exposure = [a.exposure_dose for a in acqs]
            res.resolve("dose_per_tilt", companion="per-image (companion)", note="companion")
    if acq is None:
        acq_file = res.optional("acq_order")
        if acq_file:
            t = AreTomo3TLT.from_file(str(acq_file))
            if t.n_rows != n_raw or not t.has_acq_index:
                raise ValueError(f"--acq-order {acq_file}: need {n_raw} rows with an acquisition index column")
            acq = list(t.acq_indices)
            if t.has_dose:
                exposure = list(t.doses)
    if acq is not None and exposure is None:
        d = res.optional("dose_per_tilt")
        if d is not None:
            exposure = [float(d)] * n_raw
    if acq is not None and exposure is not None:
        tlt = AreTomo3TLT.from_aln(aln, acq_index_1b=acq, dose=exposure)
    else:
        tlt = AreTomo3TLT.from_aln(aln)
        sr.warnings.append(
            "no acquisition order / exposure (companion, --acq-order, --dose-per-tilt): _TLT.txt written with the tilt column only"
        )

    # _CTF.txt when every image carries CTF metadata
    ctf_file = None
    ctfs = [im.ctf_metadata for im in images]
    if all(c is not None and c.defocus_u is not None for c in ctfs) and not res.optional("no_ctf", absent=False):
        df_hand = res.optional("defocus_hand", companion=(comp_ts.defocus_hand if comp_ts else None))
        rows = [
            cets_ctf.to_aretomo3_row(c, i + 1, df_hand=None if df_hand is None else int(df_hand))
            for i, c in enumerate(ctfs)
        ]
        ctf_file = AreTomo3CTF(rows=rows)
        if any(c.phase_shift is None for c in ctfs):
            sr.warnings.append("some images have no phase_shift: 0 written")
    elif any(c is not None for c in ctfs):
        sr.warnings.append("CTF metadata is incomplete across the tilt images: no _CTF.txt written")

    # write
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: Dict[str, Path] = {"aln": out_dir / f"{stem}.aln", "tlt": out_dir / f"{stem}_TLT.txt"}
    if ctf_file is not None:
        outputs["ctf"] = out_dir / f"{stem}_CTF.txt"
    for p in outputs.values():
        if p.exists() and not overwrite:
            raise FileExistsError(f"{p} exists (use --overwrite)")
    aln.to_file(str(outputs["aln"]))
    tlt.to_file(str(outputs["tlt"]))
    if ctf_file is not None:
        ctf_file.to_file(str(outputs["ctf"]))

    stack_src = _resolve_doc_path(ts.path, doc_dir)
    if res.optional("link_stack", absent=True) and stack_src is not None and stack_src.exists():
        link = out_dir / f"{stem}.mrc"
        if link.is_symlink() or (link.exists() and overwrite):
            link.unlink()
        if not link.exists():
            os.symlink(stack_src.resolve(), link)
        outputs["stack"] = link
    elif stack_src is not None:
        sr.warnings.append(
            f"tilt stack {ts.path} not found next to the document; place it at {out_dir / (stem + '.mrc')}"
        )

    # gates: re-read what was written
    from cryoet_alignment.io.aretomo3 import AreTomo3ALN

    reread = AreTomo3ALN.from_file(str(outputs["aln"]))
    sr.gates.append(
        Gate(
            "aln_invariants",
            reread.n_raw == n_raw and len(reread.GlobalAlignments) == len(aligned),
            value=(reread.n_raw, len(reread.GlobalAlignments)),
            expected=(n_raw, len(aligned)),
        )
    )
    sr.gates.append(Gate("tlt_rows", AreTomo3TLT.from_file(str(outputs["tlt"])).n_rows == n_raw, expected=n_raw))
    sr.hints.append(
        aretomo3_hint(
            out_dir,
            stem,
            pixel_size_a=pix,
            voltage_kv=res.optional("voltage", companion=(comp_ts.voltage_kv if comp_ts else None)),
            cs_mm=res.optional("cs", companion=(comp_ts.cs_mm if comp_ts else None)),
            amplitude_contrast=res.optional(
                "amp_contrast", companion=(comp_ts.amplitude_contrast if comp_ts else None)
            ),
            vol_z_px=int(vol_z),
            has_ctf=ctf_file is not None,
        )
    )
    sr.outputs.update({k: str(v) for k, v in outputs.items()})
    sr.provenance = res.provenance()
    return outputs
