"""A synthetic AreTomo3 run directory built from real text fixtures (24jul16a Position_16_3: .aln with 16
patches, mdoc, _CTF.txt) plus zero-filled MRC headers sized to a 512x512 stack (the .aln RawSize is
rewritten accordingly; TX/TY are left as they are — they are just numbers to round-trip)."""

import re
import shutil
from pathlib import Path

import mrcfile
import numpy as np
import pytest

DATA = Path(__file__).parent / "data"
STEM = "Position_16_3"
NX, NY, NZ = 512, 512, 31
PIX = 1.54
BIN = 4


def write_mrc(path: Path, shape_zyx, voxel: float) -> None:
    with mrcfile.new(str(path), overwrite=True) as m:
        m.set_data(np.zeros(shape_zyx, dtype=np.float32))
        m.voxel_size = voxel


def make_run(root: Path, *, darks=(), with_vol=True, xzy=False, mdoc=True, ctf=True, exposures=None) -> Path:
    """Write a run directory under ``root`` and return its path.

    ``darks``: raw sections to turn into DarkFrame lines (their global row is removed; SEC of the
    others is preserved, as AreTomo3 does). ``exposures``: per raw section, replaces the mdoc dose
    with a ``_TLT.txt`` (acquisition order from the mdoc)."""
    run = root / "aretomo3"
    run.mkdir(parents=True, exist_ok=True)
    text = (DATA / f"{STEM}.aln").read_text()
    text = re.sub(r"# RawSize = \d+ \d+ \d+", f"# RawSize = {NX} {NY} {NZ}", text)
    if darks:
        lines = text.splitlines()
        out = []
        globals_seen = 0
        for ln in lines:
            if ln.startswith("# AlphaOffset"):
                for z in darks:
                    tilt = -45.01 + 3.0 * z
                    out.append(f"# DarkFrame =  {z:4d} {z + 1:4d} {tilt:8.2f}")
            if re.match(r"^\s+\d+\s+[-\d.]+\s+[\d.]+\s", ln) and not ln.startswith("# "):
                sec = int(ln.split()[0])
                if sec - 1 in darks and globals_seen < NZ:
                    globals_seen += 1
                    continue
                globals_seen += 1
            out.append(ln)
        text = "\n".join(out) + "\n"
        # drop the local rows of the dark tilts (sec_idx over the dark-removed list): simplest is to drop
        # all local rows and set NumPatches 0 for dark-frame fixtures
        head, _, _ = text.partition("# Local Alignment\n")
        text = head.replace("# NumPatches = 16", "# NumPatches = 0")
    (run / f"{STEM}.aln").write_text(text)
    write_mrc(run / f"{STEM}.mrc", (NZ, NY, NX), PIX)
    if with_vol:
        if xzy:
            write_mrc(run / f"{STEM}_Vol.mrc", (NY // BIN, 63, NX // BIN), PIX * BIN)  # nz = y extent
        else:
            write_mrc(run / f"{STEM}_Vol.mrc", (63, NY // BIN, NX // BIN), PIX * BIN)
    if mdoc:
        shutil.copy(DATA / f"{STEM}.mdoc", run / f"{STEM}.mdoc")
    if ctf:
        shutil.copy(DATA / f"{STEM}_CTF.txt", run / f"{STEM}_CTF.txt")
    if exposures is not None:
        from cets_aretomo3.discover import discover

        d = discover(run / f"{STEM}.aln")
        acq = d.get("acq_index_1b")
        stage = d.get("stage_tilt_deg")
        (run / f"{STEM}_TLT.txt").write_text(
            "".join(f"{stage[i]:8.2f}  {acq[i]:4d}  {exposures[i]:8.2f}\n" for i in range(NZ)),
        )
    return run


@pytest.fixture
def run_dir(tmp_path) -> Path:
    return make_run(tmp_path)
