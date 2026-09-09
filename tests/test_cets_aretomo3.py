import json
from pathlib import Path

import pytest
from click.testing import CliRunner
from conftest import BIN, NX, NZ, PIX, STEM, make_run
from cryoet_alignment.io.aretomo3 import AreTomo3ALN, AreTomo3CTF, AreTomo3TLT
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.entities import load_dataset

from cets_aretomo3.cli import main


def _run(args, expect_ok=True):
    result = CliRunner().invoke(main, args, catch_exceptions=False)
    if expect_ok:
        assert result.exit_code == 0, result.stdout + result.stderr
    return result


def _to_cets(run: Path, out: Path, *extra):
    return _run(["to-cets", str(run), "-o", str(out), "--drop-locals", *extra])


def test_to_cets_document_and_companion(run_dir, tmp_path):
    out = tmp_path / "out" / "run.cets.json"
    r = _to_cets(run_dir, out, "--flip-vol", "1")
    assert "pix = 1.54  [discovered]" in r.stdout
    assert "WARNING: dose_convention defaulted to 'exclusive'" in r.stderr
    assert "WARNING: paths defaulted to 'relative'" in r.stderr
    ds = load_dataset(out)
    region = ds.regions[0]
    ts = region.tilt_series[0]
    assert ts.id == STEM and len(ts.images) == NZ
    assert ts.images[0].nominal_tilt_angle == pytest.approx(-45.01)
    assert ts.images[0].ctf_metadata.defocus_u == pytest.approx(21717.31)
    assert ts.images[0].ctf_metadata.defocus_handedness is None  # dfHand is companion-only
    assert ts.images[0].width == NX and ts.images[0].coordinate_transformations[0].sequence[1].scale == [PIX, PIX]
    # accumulated dose = exclusive pre-exposure from the mdoc order/exposure
    comp = Companion.load_for(out)
    imgs = comp.tilt_series[STEM].images
    acq = {i: imgs[f"{STEM}_{i}"].acquisition_index_1b for i in range(NZ)}
    exp = {i: imgs[f"{STEM}_{i}"].exposure_dose for i in range(NZ)}
    for i in range(NZ):
        before = sum(exp[j] for j in range(NZ) if acq[j] < acq[i])
        assert ts.images[i].accumulated_dose == pytest.approx(before)
    assert comp.tilt_series[STEM].defocus_hand == -1
    assert comp.tilt_series[STEM].voltage_kv == 300.0
    # tomograms: the bin-1 reference box and the declared _Vol.mrc
    ids = [t.id for t in region.tomograms]
    assert ids == [f"{STEM}_volume", f"{STEM}_tomo"]
    assert (region.tomograms[0].width, region.tomograms[0].depth) == (NX, 63 * BIN)
    assert region.tomograms[1].coordinate_transformations[0].sequence[1].scale == [PIX * BIN] * 3
    a = comp.alignments[0]
    assert a.reference_tomogram_id == f"{STEM}_volume" and a.tomogram_ids == ids
    assert a.thickness_px == 1500 and a.alignment_type == "GLOBAL" and "local alignment" in a.dropped[0]
    assert comp.tomograms[f"{STEM}_tomo"].flip_vol == 1
    assert region.movie_stack_collection.movie_stacks[0].stacks[0].path.endswith(".eer")
    report = json.loads((out.parent / "run.cets.report.json").read_text())
    assert report["ok"] and report["series"][0]["gates"][0]["name"] == "rows"


def test_locals_refused_without_flag(run_dir, tmp_path):
    r = _run(["to-cets", str(run_dir), "-o", str(tmp_path / "x.cets.json")], expect_ok=False)
    assert r.exit_code == 1 and "local alignment" in r.stderr


def test_flip_vol_never_inferred(tmp_path):
    run = make_run(tmp_path)
    out = tmp_path / "o" / "a.cets.json"
    r = _to_cets(run, out)
    assert "FlipVol not declared" in r.stderr
    assert [t.id for t in load_dataset(out).regions[0].tomograms] == [f"{STEM}_volume"]
    run2 = make_run(tmp_path / "xzy", xzy=True)
    r = _to_cets(run2, tmp_path / "o2" / "b.cets.json", "--flip-vol", "1")
    assert "XZY" in r.stderr and "-FlipVol 0" in r.stderr


def test_roundtrip_aln_and_tlt_identical(run_dir, tmp_path):
    out = tmp_path / "out" / "run.cets.json"
    _to_cets(run_dir, out, "--flip-vol", "1", "--voltage", "300", "--cs", "2.7", "--amp-contrast", "0.07")
    back = tmp_path / "back"
    r = _run(["from-cets", str(out), "-o", str(back)])
    assert "-Cmd 2" in r.stdout and "-CorrCTF 1" in r.stdout and "-VolZ 252" in r.stdout
    src = AreTomo3ALN.from_file(run_dir / f"{STEM}.aln")
    got = AreTomo3ALN.from_file(back / f"{STEM}.aln")
    assert [str(g) for g in got.GlobalAlignments] == [str(g) for g in src.GlobalAlignments]
    assert got.NumPatches == 0 and got.Thickness == src.Thickness == 1500
    tlt = AreTomo3TLT.from_file(back / f"{STEM}_TLT.txt")
    comp = Companion.load_for(out)
    assert tlt.acq_indices == [comp.tilt_series[STEM].images[f"{STEM}_{i}"].acquisition_index_1b for i in range(NZ)]
    assert tlt.doses == pytest.approx([comp.tilt_series[STEM].images[f"{STEM}_{i}"].exposure_dose for i in range(NZ)])
    ctf = AreTomo3CTF.from_file(back / f"{STEM}_CTF.txt")
    src_ctf = AreTomo3CTF.from_file(run_dir / f"{STEM}_CTF.txt")
    for a, b in zip(src_ctf.rows, ctf.rows, strict=True):
        # the profile wraps the astigmatism angle into [0, 180) (an identity of the cos(2(theta-angle)) field)
        assert (a.df_max_a, a.df_min_a, a.azimuth_deg % 180.0) == pytest.approx((b.df_max_a, b.df_min_a, b.azimuth_deg), abs=1e-6)
        assert a.phase_rad == pytest.approx(b.phase_rad, abs=1e-9) and b.df_hand == -1  # from the companion
    assert (back / f"{STEM}.mrc").is_symlink()


def test_dark_frames_roundtrip(tmp_path):
    run = make_run(tmp_path, darks=(2, 29), mdoc=False, ctf=False)
    src = AreTomo3ALN.from_file(run / f"{STEM}.aln")
    assert src.dark_indices() == [2, 29] and len(src.GlobalAlignments) == NZ - 2
    out = tmp_path / "o" / "d.cets.json"
    _run(["to-cets", str(run), "-o", str(out), "--tomo-size", "300"])
    ds = load_dataset(out)
    region = ds.regions[0]
    assert len(region.tilt_series[0].images) == NZ
    assert [pa.tilt_image_id for pa in region.alignments[0].projection_alignments] == [f"{STEM}_{z}" for z in range(NZ) if z not in (2, 29)]
    back = tmp_path / "back"
    r = _run(["from-cets", str(out), "-o", str(back)])
    assert "tilt column only" in r.stderr  # no order/exposure without mdoc/_TLT.txt
    got = AreTomo3ALN.from_file(back / f"{STEM}.aln")
    assert [str(d) for d in got.DarkFrames] == [str(d) for d in src.DarkFrames]
    assert [g.sec for g in got.GlobalAlignments] == [g.sec for g in src.GlobalAlignments]
    assert AreTomo3TLT.from_file(back / f"{STEM}_TLT.txt").tilts == pytest.approx(src.raw_tilts())


def test_exposure_preserved_only_with_companion(tmp_path):
    exposures = [2.0, 3.0, 100.0] + [1.0] * (NZ - 3)
    run = make_run(tmp_path, ctf=False, exposures=exposures)
    out = tmp_path / "o" / "e.cets.json"
    _to_cets(run, out, "--flip-vol", "1")
    src_tlt = AreTomo3TLT.from_file(run / f"{STEM}_TLT.txt")
    back = tmp_path / "back"
    _run(["from-cets", str(out), "-o", str(back)])
    got = AreTomo3TLT.from_file(back / f"{STEM}_TLT.txt")
    assert got.doses == pytest.approx(src_tlt.doses) and got.acq_indices == src_tlt.acq_indices
    # without the companion the exposures are gone: 1-column file + warning, or --dose-per-tilt + --acq-order
    Companion.path_for(out).unlink()
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "back2")], expect_ok=False)
    assert "select one with --tomogram" in r.stderr  # the companion also carried the reference-tomogram binding
    ref = ["--tomogram", f"{STEM}_volume"]
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "back2"), *ref])
    assert "tilt column only" in r.stderr
    assert not AreTomo3TLT.from_file(tmp_path / "back2" / f"{STEM}_TLT.txt").has_acq_index
    _run(["from-cets", str(out), "-o", str(tmp_path / "back3"), *ref, "--acq-order", str(run / f"{STEM}_TLT.txt")])
    assert AreTomo3TLT.from_file(tmp_path / "back3" / f"{STEM}_TLT.txt").doses == pytest.approx(src_tlt.doses)


def test_x_rotation_refused_unless_dropped(run_dir, tmp_path):
    out = tmp_path / "o" / "x.cets.json"
    _to_cets(run_dir, out, "--flip-vol", "1")
    doc = json.loads(out.read_text())
    import numpy as np
    from cryoet_alignment.io.cets.rotation import tilt_matrix

    for pa in doc["regions"][0]["alignments"][0]["projection_alignments"]:
        tilt = float(np.degrees(np.arcsin(-pa["sequence"][0]["affine"][2][0])))
        pa["sequence"][0]["affine"] = tilt_matrix(tilt, 0.35).tolist()
    out.write_text(json.dumps(doc))
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "b")], expect_ok=False)
    assert r.exit_code == 1 and "X rotation" in r.stderr
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "b"), "--drop-x-rotation"])
    assert "dropped: volume X rotation up to 0.35" in r.stderr


def test_config_and_selection(run_dir, tmp_path):
    cfg = tmp_path / "cets.yaml"
    cfg.write_text("cets:\n  dose_convention: exclusive\n  paths: absolute\ncets-aretomo3:\n  to-cets: {flip_vol: 1}\n")
    out = tmp_path / "o" / "c.cets.json"
    r = _to_cets(run_dir, out, "--config", str(cfg))
    assert "defaulted" not in r.stderr
    assert "paths = 'absolute'  [config]" in r.stdout and "flip_vol = 1  [config]" in r.stdout
    assert Path(load_dataset(out).regions[0].tilt_series[0].path).is_absolute()
    (tmp_path / "bad.yaml").write_text("cets:\n  pixl: 1\n")
    r = _run(["to-cets", str(run_dir), "-o", str(tmp_path / "z.cets.json"), "--config", str(tmp_path / "bad.yaml")], expect_ok=False)
    assert "unknown option" in (r.stderr + r.stdout)
    # two alignments in one region -> explicit selection
    doc = json.loads(out.read_text())
    reg = doc["regions"][0]
    dup = json.loads(json.dumps(reg["alignments"][0]))
    for pa in dup["projection_alignments"]:
        pa["id"] = pa["id"].replace("_aretomo3_align_", "_other_align_")
    reg["alignments"].append(dup)
    out.write_text(json.dumps(doc))
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "s")], expect_ok=False)
    assert "select one with --alignment" in r.stderr
    _run(["from-cets", str(out), "-o", str(tmp_path / "s"), "--alignment", "other", "--tomogram", f"{STEM}_volume"])
