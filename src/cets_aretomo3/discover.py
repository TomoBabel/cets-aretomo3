"""Discover a series' metadata from the files AreTomo3 leaves next to its ``.aln`` (header-only reads) and
the mdoc when present. Every discovered value is reported with its source; explicit values always win.

Never consulted: ``AreTomo3_Session.json``. Never inferred: FlipVol 1 vs 2 (both write XYZ volumes,
``CAreTomoMain.cpp:893``), defocus hand, voltage without an mdoc.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import mrcfile
import numpy as np
from cryoet_alignment.io.aretomo3 import AreTomo3ALN, AreTomo3CTF, AreTomo3TLT

_TILT_AXIS_RE = re.compile(r"TiltAxisAngle\s*=\s*([-+]?\d+(?:\.\d+)?)")
STACK_EXTS = (".mrc", ".st", ".mrcs")


@dataclass
class Found:
    value: Any
    source: str


@dataclass
class AreTomo3Run:
    """Everything discovered for one series (values are ``Found(value, source)``)."""

    stem: str
    aln_path: Path
    aln: AreTomo3ALN
    stack_path: Optional[Path] = None
    tlt_path: Optional[Path] = None
    ctf_path: Optional[Path] = None
    mdoc_path: Optional[Path] = None
    vol_path: Optional[Path] = None
    found: Dict[str, Found] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def get(self, key: str, default=None):
        f = self.found.get(key)
        return default if f is None else f.value

    def source(self, key: str) -> Optional[str]:
        f = self.found.get(key)
        return None if f is None else f.source

    def add(self, key: str, value: Any, source: str, *, override: bool = False) -> None:
        if key not in self.found or override:
            self.found[key] = Found(value, source)


def mrc_header(path: Path) -> dict:
    with mrcfile.open(str(path), header_only=True, permissive=True) as m:
        h = m.header
        vs = m.voxel_size
        return {
            "nx": int(h.nx),
            "ny": int(h.ny),
            "nz": int(h.nz),
            "mode": int(h.mode),
            "voxel": (float(vs.x), float(vs.y), float(vs.z)),
        }


def _find_stack(stem: str, dirs: List[Path]) -> Optional[Path]:
    for d in dirs:
        for ext in STACK_EXTS:
            p = d / f"{stem}{ext}"
            if p.exists():
                return p
    return None


def _mdoc_sections(mdoc_path: Path) -> Tuple[Any, List[dict]]:
    from mdocfile.data_models import Mdoc

    m = Mdoc.from_file(str(mdoc_path))
    out = []
    for s in m.section_data:
        if s.ZValue is None:
            continue
        d = {k: v for k, v in s.model_dump().items() if v is not None}
        d["ZValue"] = int(s.ZValue)
        out.append(d)
    return m, out


def _time_key(entry: dict):
    from datetime import datetime

    text = entry.get("DateTime")
    if text:
        s = str(text).strip()
        for fmt in ("%d-%b-%Y %H:%M:%S", "%d-%b-%y %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%d-%b-%Y  %H:%M:%S"):
            try:
                return (0, datetime.strptime(s, fmt).timestamp(), int(entry["ZValue"]))  # noqa: DTZ007
            except ValueError:
                continue
        return (1, s, int(entry["ZValue"]))
    return (2, 0.0, int(entry["ZValue"]))


def mdoc_acquisition(sections: list, raw_stage_tilts: np.ndarray, tol_deg: float = 0.5):
    """Per raw row: (1-based acquisition index, mdoc section) matched by stage angle after removing the mean
    offset between the two angle lists."""
    if len(sections) != len(raw_stage_tilts):
        raise ValueError(f"mdoc has {len(sections)} sections for {len(raw_stage_tilts)} raw tilts")
    ordered = sorted(sections, key=_time_key)
    mtilt = np.array([float(e["TiltAngle"]) for e in ordered], dtype=np.float64)
    raw = np.asarray(raw_stage_tilts, dtype=np.float64)
    mtilt = mtilt + float(raw.mean() - mtilt.mean())
    used = np.zeros(len(ordered), dtype=bool)
    acq, rows = [], []
    for tilt in raw:
        diff = np.abs(mtilt - tilt) + used * 1e6
        j = int(np.argmin(diff))
        if diff[j] > tol_deg:
            raise ValueError(f"mdoc: no section within {tol_deg} deg of raw tilt {tilt:.2f}")
        used[j] = True
        acq.append(j + 1)
        rows.append(ordered[j])
    return acq, rows


def _sub_frame_name(value) -> str:
    return Path(str(value).replace("\\", "/")).name


def discover(
    aln_path: Path,
    *,
    stack: Optional[Path] = None,
    tilt_stack_dir: Optional[Path] = None,
    mdoc: Optional[Path] = None,
    mdoc_dir: Optional[Path] = None,
    tlt: Optional[Path] = None,
    ctf: Optional[Path] = None,
    no_ctf: bool = False,
    adjacent: bool = True,
) -> AreTomo3Run:
    aln_path = Path(aln_path)
    stem = aln_path.stem
    aln = AreTomo3ALN.from_file(str(aln_path))
    run = AreTomo3Run(stem=stem, aln_path=aln_path, aln=aln)
    n_raw = aln.n_raw
    alpha = float(aln.AlphaOffset or 0.0)

    # stage angles per raw row derived from the .aln: TILT - AlphaOffset on global rows, DarkFrame angle on darks
    stage = np.array(aln.raw_tilts(), dtype=np.float64)
    for g in aln.GlobalAlignments:
        stage[g.sec - 1] = g.tilt - alpha
    run.add("stage_tilt_deg", stage.tolist(), f"{aln_path.name}#TILT-AlphaOffset")
    run.add("image_dims_px", (int(aln.RawSize[0]), int(aln.RawSize[1])), f"{aln_path.name}#RawSize")
    run.add("n_raw", n_raw, f"{aln_path.name}#RawSize")
    run.add("alpha_offset_deg", alpha, f"{aln_path.name}#AlphaOffset")
    run.add("beta_offset_deg", float(aln.BetaOffset or 0.0), f"{aln_path.name}#BetaOffset")
    if aln.Thickness:
        run.add("thickness_px", int(aln.Thickness), f"{aln_path.name}#Thickness")
    if aln.GlobalAlignments:
        run.add("tilt_axis_deg", float(np.median([g.rot for g in aln.GlobalAlignments])), f"{aln_path.name}#median ROT")

    # tilt stack: dims + pixel size from the header
    stack_path = Path(stack) if stack else None
    if stack_path is None:
        dirs = [Path(tilt_stack_dir)] if tilt_stack_dir else []
        if adjacent:
            dirs.append(aln_path.parent)
        stack_path = _find_stack(stem, dirs) if dirs else None
    if stack_path is not None and stack_path.exists():
        h = mrc_header(stack_path)
        run.stack_path = stack_path
        if (h["nx"], h["ny"]) != tuple(run.get("image_dims_px")):
            run.warnings.append(
                f"{stack_path.name}: header {h['nx']}x{h['ny']} differs from .aln RawSize {run.get('image_dims_px')}"
            )
        if h["nz"] != n_raw:
            run.warnings.append(f"{stack_path.name}: {h['nz']} sections for RawSize z = {n_raw}")
        run.add("image_dims_px", (h["nx"], h["ny"]), f"{stack_path.name}#header", override=True)
        if h["voxel"][0] > 0:
            run.add("pix", round(h["voxel"][0], 6), f"{stack_path.name}#header")

    # _TLT.txt: stage angles, acquisition order, per-image dose
    tlt_path = Path(tlt) if tlt else (aln_path.with_name(f"{stem}_TLT.txt") if adjacent else None)
    if tlt_path is not None and tlt_path.exists():
        t = AreTomo3TLT.from_file(str(tlt_path))
        if t.n_rows != n_raw:
            run.warnings.append(f"{tlt_path.name}: {t.n_rows} rows for {n_raw} raw sections - ignored")
        else:
            run.tlt_path = tlt_path
            dev = float(np.abs(np.array(t.tilts) - stage).max())
            if dev > 0.5:
                run.warnings.append(
                    f"{tlt_path.name}: tilt column deviates from .aln stage angles by up to {dev:.2f} deg"
                )
            run.add("stage_tilt_deg", [float(v) for v in t.tilts], tlt_path.name, override=True)
            if t.has_acq_index:
                run.add("acq_index_1b", list(t.acq_indices), tlt_path.name)
            if t.has_dose and sum(t.doses) > 0:
                run.add("exposure", list(t.doses), tlt_path.name)

    # mdoc: pixel size, voltage, order, dose, frame names
    mdoc_path = Path(mdoc) if mdoc else None
    if mdoc_path is None:
        cands = ([Path(mdoc_dir) / f"{stem}.mdoc"] if mdoc_dir else []) + (
            [aln_path.with_name(f"{stem}.mdoc")] if adjacent else []
        )
        mdoc_path = next((c for c in cands if c.exists()), None)
    if mdoc_path is not None and mdoc_path.exists():
        m, sections = _mdoc_sections(mdoc_path)
        run.mdoc_path = mdoc_path
        g = m.global_data
        pix = g.PixelSpacing or (sections[0].get("PixelSpacing") if sections else None)
        if pix:
            run.add("pix", float(pix), f"{mdoc_path.name}#PixelSpacing")
        volt = g.Voltage or (sections[0].get("Voltage") if sections else None)
        if volt:
            run.add("voltage", float(volt), f"{mdoc_path.name}#Voltage")
        for title in m.titles or []:
            mt = _TILT_AXIS_RE.search(title)
            if mt:
                run.add("tilt_axis_nominal_deg", float(mt.group(1)), f"{mdoc_path.name}#TiltAxisAngle")
        try:
            acq, rows = mdoc_acquisition(sections, stage)
        except ValueError as e:
            run.warnings.append(str(e))
        else:
            # _TLT.txt (AreTomo3's own record of order/dose/angles) wins over the mdoc when both exist
            run.add("acq_index_1b", acq, mdoc_path.name)
            run.add(
                "stage_tilt_deg", [float(r["TiltAngle"]) for r in rows], mdoc_path.name, override=run.tlt_path is None
            )
            doses = [float(r.get("ExposureDose", 0.0) or 0.0) for r in rows]
            if sum(doses) > 0:
                run.add("exposure", doses, f"{mdoc_path.name}#ExposureDose")
            if all(r.get("SubFramePath") for r in rows):
                run.add(
                    "frame_names", [_sub_frame_name(r["SubFramePath"]) for r in rows], f"{mdoc_path.name}#SubFramePath"
                )

    # CTF
    if not no_ctf:
        ctf_path = Path(ctf) if ctf else (aln_path.with_name(f"{stem}_CTF.txt") if adjacent else None)
        if ctf_path is not None and ctf_path.exists():
            c = AreTomo3CTF.from_file(str(ctf_path))
            if c.n_rows != n_raw:
                run.warnings.append(f"{ctf_path.name}: {c.n_rows} rows for {n_raw} raw sections - ignored")
            else:
                run.ctf_path = ctf_path
                run.add("ctf", c, ctf_path.name)

    # tomogram: dims + bin from <stem>_Vol.mrc (FlipVol 1 vs 2 is NOT inferable from the header)
    vol = aln_path.with_name(f"{stem}_Vol.mrc")
    if adjacent and vol.exists():
        try:
            h = mrc_header(vol)
        except Exception as e:  # noqa: BLE001
            run.warnings.append(f"{vol.name}: unreadable header ({e})")
        else:
            run.vol_path = vol
            rx, ry = run.get("image_dims_px")
            b = rx / h["nx"] if h["nx"] else 0.0
            if b > 0 and abs(h["ny"] * b - ry) <= b:
                run.add("vol_layout", "xyz", f"{vol.name}#header")  # FlipVol 1 or 2
                run.add("tomo_dims_px", (h["nx"], h["ny"], h["nz"]), f"{vol.name}#header")
            elif b > 0 and abs(h["nz"] * b - ry) <= b:
                run.add("vol_layout", "xzy", f"{vol.name}#header")  # FlipVol 0
                run.add("tomo_dims_px", (h["nx"], h["nz"], h["ny"]), f"{vol.name}#header (xzy)")
            else:
                run.warnings.append(
                    f"{vol.name}: header {h['nx']}x{h['ny']}x{h['nz']} does not match RawSize {rx}x{ry} at any bin"
                )
            if b > 0:
                run.add("bin", b, f"{vol.name}#header / RawSize")
                run.add("vol_voxel_header_a", round(h["voxel"][0], 6), f"{vol.name}#header")  # float32 header
    return run
