"""``cets-aretomo3`` command line: ``to-cets`` (AreTomo3 run -> CETS) and ``from-cets`` (CETS -> -Cmd 2 inputs)."""

import sys
from pathlib import Path

import click
from cryoet_alignment.io.cets.cli_support import (
    Report,
    SeriesReport,
    common_options,
    echo,
    expand_sources,
    finish,
    load_config,
    make_resolver,
    parse_alignment_selector,
    print_series,
    selection_options,
)
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.entities import dataset_entity, dump_json, load_dataset, validate_document

from cets_aretomo3 import __version__
from cets_aretomo3.discover import discover
from cets_aretomo3.from_cets import cets_to_aretomo3
from cets_aretomo3.to_cets import aretomo3_to_cets

PACKAGE = "cets-aretomo3"
TO_CETS_OPTIONS = {
    "pix",
    "tomo_size",
    "flip_vol",
    "mdoc_dir",
    "tilt_stack_dir",
    "dose_convention",
    "dose_per_tilt",
    "drop_locals",
    "no_ctf",
    "paths",
    "voltage",
    "cs",
    "amp_contrast",
}
FROM_CETS_OPTIONS = {
    "tomo_size",
    "drop_x_rotation",
    "dose_per_tilt",
    "acq_order",
    "defocus_hand",
    "link_stack",
    "no_ctf",
    "voltage",
    "cs",
    "amp_contrast",
}


@click.group()
@click.version_option(__version__)
def main():
    """AreTomo3 <-> CETS (rigid profile cets-rigid/0.1)."""


@main.command("to-cets")
@click.argument("sources", nargs=-1, required=True)
@click.option(
    "-o", "--output", "output", required=True, type=click.Path(dir_okay=False), help="CETS dataset JSON to write."
)
@click.option("--name", default=None, help="Dataset name (default: the output stem).")
@click.option("--pix", type=float, default=None, help="Tilt-image pixel size (Å/px).")
@click.option(
    "--tomo-size", "tomo_size", type=float, default=None, help="Reconstruction thickness Z in unbinned pixels."
)
@click.option(
    "--flip-vol",
    "flip_vol",
    type=click.Choice(["0", "1", "2"]),
    default=None,
    help="The -FlipVol the run used (needed to describe <stem>_Vol.mrc; never inferred).",
)
@click.option(
    "--mdoc-dir",
    "mdoc_dir",
    type=click.Path(file_okay=False, exists=True),
    default=None,
    help="Directory holding <stem>.mdoc (default: next to the .aln).",
)
@click.option(
    "--tilt-stack-dir",
    "tilt_stack_dir",
    type=click.Path(file_okay=False, exists=True),
    default=None,
    help="Directory holding <stem>.mrc/.st/.mrcs (default: next to the .aln).",
)
@click.option(
    "--dose-convention",
    "dose_convention",
    type=click.Choice(["exclusive", "inclusive"]),
    default=None,
    help="Whether accumulated_dose excludes the image's own exposure [exclusive, warned].",
)
@click.option("--dose-per-tilt", "dose_per_tilt", type=float, default=None, help="Constant per-image exposure (e/Å²).")
@click.option(
    "--drop-locals", "drop_locals", is_flag=True, default=None, help="Keep the global rows of a local alignment."
)
@click.option("--no-ctf", "no_ctf", is_flag=True, default=None, help="Ignore <stem>_CTF.txt.")
@click.option(
    "--paths",
    type=click.Choice(["relative", "absolute"]),
    default=None,
    help="How file paths are written into the document [relative, warned].",
)
@click.option("--voltage", type=float, default=None, help="kV (companion only).")
@click.option("--cs", type=float, default=None, help="mm (companion only).")
@click.option("--amp-contrast", "amp_contrast", type=float, default=None, help="Amplitude contrast (companion only).")
@common_options
def to_cets(sources, output, name, config_path, overwrite, fail_fast, **cli):
    """Convert AreTomo3 runs (.aln files or directories) to one CETS dataset JSON + companion + report."""
    out = Path(output)
    if out.exists() and not overwrite:
        raise click.ClickException(f"{out} exists (use --overwrite)")
    out.parent.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, TO_CETS_OPTIONS)
    flags = {k: v for k, v in cli.items() if v is not None}
    if "flip_vol" in flags:
        flags["flip_vol"] = int(flags["flip_vol"])
    report = Report(PACKAGE, "to-cets")
    companion = Companion(generator=f"{PACKAGE} {__version__}")
    regions = []
    for aln_path in expand_sources(sources, ["*.aln"]):
        sr = SeriesReport(aln_path.stem)
        report.series.append(sr)
        try:
            run = discover(
                aln_path,
                mdoc_dir=Path(flags["mdoc_dir"]) if "mdoc_dir" in flags else None,
                tilt_stack_dir=Path(flags["tilt_stack_dir"]) if "tilt_stack_dir" in flags else None,
                no_ctf=bool(flags.get("no_ctf")),
            )
            sr.warnings.extend(run.warnings)
            res = make_resolver(PACKAGE, "to-cets", flags, config, run.stem, sr)
            result = aretomo3_to_cets(run, res, sr, out_dir=out.parent)
        except Exception as e:  # noqa: BLE001 - reported per series
            sr.error = str(e)
            print_series(sr)
            if fail_fast:
                break
            continue
        regions.append(result.region)
        companion.tilt_series[run.stem] = result.tilt_series_companion
        companion.alignments.append(result.alignment_companion)
        companion.tomograms.update(result.tomogram_companions)
        print_series(sr)
    if regions:
        ds = dataset_entity(name or out.name.split(".")[0], regions)
        validate_document(ds)
        dump_json(ds, out)
        companion.dump(Companion.path_for(out))
        echo(f"wrote {out} ({len(regions)} region(s)) + {Companion.path_for(out).name}")
    finish(report, out.with_name(out.name.split(".")[0] + ".cets.report.json"))


@main.command("from-cets")
@click.argument("document", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o", "--output", "output", required=True, type=click.Path(file_okay=False), help="AreTomo3 -Cmd 2 input directory."
)
@selection_options
@click.option(
    "--tomo-size",
    "tomo_size",
    type=int,
    default=None,
    help="Thickness Z in unbinned pixels (default: the reference tomogram).",
)
@click.option(
    "--drop-x-rotation",
    "drop_x_rotation",
    is_flag=True,
    default=None,
    help="Write alignments that carry volume X rotations, dropping them (an .aln cannot hold them).",
)
@click.option("--dose-per-tilt", "dose_per_tilt", type=float, default=None, help="Constant per-image exposure (e/Å²).")
@click.option(
    "--acq-order",
    "acq_order",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="A _TLT.txt-style file with the acquisition index (and dose) per raw section.",
)
@click.option(
    "--defocus-hand",
    "defocus_hand",
    type=click.Choice(["-1", "1"]),
    default=None,
    help="Defocus handedness column of _CTF.txt (default: companion; omitted when unknown).",
)
@click.option(
    "--link-stack/--no-link-stack",
    "link_stack",
    default=None,
    help="Symlink the tilt stack into the output directory [on].",
)
@click.option("--no-ctf", "no_ctf", is_flag=True, default=None, help="Do not write <stem>_CTF.txt.")
@click.option("--voltage", type=float, default=None, help="kV for the -Cmd 2 hint (default: companion).")
@click.option("--cs", type=float, default=None, help="mm for the -Cmd 2 hint (default: companion).")
@click.option("--amp-contrast", "amp_contrast", type=float, default=None, help="Amplitude contrast for the hint.")
@common_options
def from_cets(document, output, regions, alignment, tomogram, config_path, overwrite, fail_fast, **cli):
    """Write AreTomo3 -Cmd 2 inputs (.aln, _TLT.txt, _CTF.txt) for the regions of a CETS dataset."""
    doc = Path(document)
    ds = load_dataset(doc)
    companion = Companion.load_for(doc)
    out_dir = Path(output)
    config = load_config(config_path, FROM_CETS_OPTIONS)
    flags = {k: v for k, v in cli.items() if v is not None}
    if "defocus_hand" in flags:
        flags["defocus_hand"] = int(flags["defocus_hand"])
    report = Report(PACKAGE, "from-cets")
    wanted = set(regions)
    for region in ds.regions:
        if wanted and region.id not in wanted:
            continue
        sr = SeriesReport(region.id)
        report.series.append(sr)
        try:
            res = make_resolver(PACKAGE, "from-cets", flags, config, region.id, sr)
            cets_to_aretomo3(
                region,
                res,
                sr,
                out_dir=out_dir,
                doc_dir=doc.parent,
                companion=companion,
                alignment_selector=parse_alignment_selector(alignment),
                tomogram_selector=tomogram,
                overwrite=overwrite,
            )
        except Exception as e:  # noqa: BLE001
            sr.error = str(e)
        print_series(sr)
        if sr.error and fail_fast:
            break
    out_dir.mkdir(parents=True, exist_ok=True)
    finish(report, out_dir / "cets_aretomo3.report.json")


if __name__ == "__main__":
    sys.exit(main())
