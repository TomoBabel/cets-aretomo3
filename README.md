# cets-aretomo3

Convert AreTomo3 tilt-series alignments to the CETS cryo-ET standard and back.

The package implements the rigid profile `cets-rigid/0.1` of CETS, documented in
[`cryoet-alignment/docs/cets.md`](https://github.com/uermel/cryoet-alignment/blob/uermel/cets/docs/cets.md).
It handles global (rigid) alignments only. AreTomo3 patch-tracking locals are refused unless you ask
to drop them.

- `cets-aretomo3 to-cets` reads one or more AreTomo3 output directories and writes one CETS dataset.
- `cets-aretomo3 from-cets` writes the input files for an AreTomo3 `-Cmd 2` reconstruction from a CETS dataset.

Geometry (rotation operators, shifts, CTF values) round-trips through CETS alone. Acquisition order,
per-image exposure, `-FlipVol`, kV/Cs and the defocus hand are not part of CETS yet. They travel in a
companion manifest written next to every document, and the writer asks for them when it is missing.

## Installation

```bash
pip install git+https://github.com/TomoBabel/cets-aretomo3.git
```

The package depends on `cryoet-alignment >= 0.3.0` (the CETS codec) and on the pinned
`cets_data_model` commit listed in `pyproject.toml`. Until cryoet-alignment 0.3.0 is on PyPI, install it
from its branch first:

```bash
pip install git+https://github.com/uermel/cryoet-alignment.git@uermel/cets
```

## Expected source layout

`to-cets` takes `.aln` files or directories. In a directory every `*.aln` is one tilt series, and the
other files of that series are found next to it by stem.

```
aretomo3_out/                    # one directory, any number of series
├── TS_01.aln                    # required: AreTomo3 alignment (global rows; locals need --drop-locals)
├── TS_01.mrc                    # tilt stack AreTomo3 aligned (.mrc/.st/.mrcs): width, height, pixel size
│                                #   optional: else --pix or the mdoc PixelSpacing; --tilt-stack-dir DIR looks elsewhere
├── TS_01_TLT.txt                # stage angle, acquisition index, dose per raw section (AreTomo3 output)
│                                #   optional: wins over the mdoc when both exist
├── TS_01_CTF.txt                # per-section CTF (AreTomo3 output); optional, skipped with --no-ctf
├── TS_01_Vol.mrc                # reconstruction; optional: Z extent x bin gives the box thickness
│                                #   described as a Tomogram only with --flip-vol 1|2 (never inferred)
├── TS_01.mdoc                   # SerialEM mdoc: frame names, stage angles, order, exposure, voltage
│                                #   optional: or --mdoc-dir DIR with DIR/TS_01.mdoc
├── TS_02.aln
└── ...
```

Required per series: the `.aln`, a pixel size (stack header, mdoc or `--pix`) and the box thickness
(`_Vol.mrc`, or `--tomo-size` in unbinned pixels). Without `_TLT.txt` or mdoc the tilt series still
converts, but `accumulated_dose` stays null and the nominal angles come from the `.aln` itself
(`TILT − AlphaOffset` on aligned rows, the DarkFrame angle on excluded ones).

## AreTomo3 to CETS

```bash
cets-aretomo3 to-cets aretomo3_out/ -o cets/at3.cets.json --flip-vol 1
```

What happens per series:

1. The `.aln` is parsed with its invariants enforced (SEC is a permutation over the raw sections,
   DarkFrame lines agree with the missing rows).
2. Stage angles, acquisition order, exposure and frame names come from `_TLT.txt`, then the mdoc.
3. Every raw section becomes a `TiltImage` (`section` = raw z index). Sections without a global row
   (dark frames) get no `ProjectionAlignment`.
4. The alignment is expressed in the tilt image's and the reference volume's centred physical frames,
   with AreTomo3's half-pixel centre convention corrected on the way in.
5. The reference volume is the alignment's own bin-1 box (`RawSize` x pixel, thickness from `_Vol.mrc`
   or `--tomo-size`) as a virtual tomogram `<stem>_volume`. With `--flip-vol 1|2` the `_Vol.mrc` file is
   described as a second tomogram `<stem>_tomo` at its header voxel size.

Console output (abridged): every value is printed with where it came from.

```
== Position_16_3
   pix = 1.54  [discovered]  (Position_16_3.mrc#header)
   paths = 'relative'  [default]
   dose_convention = 'exclusive'  [default]
   stage_tilt_deg = [-45.01, -42.01, ...]  [discovered]  (Position_16_3.mdoc)
   acq_index_1b = [31, 28, 27, ...]  [discovered]  (Position_16_3.mdoc)
   exposure = [1.49, 1.57, ...]  [discovered]  (Position_16_3.mdoc#ExposureDose)
   ctf = 'Position_16_3_CTF.txt'  [discovered]  (_CTF.txt)
   tomo_size = 252.0  [discovered]  (Position_16_3_Vol.mrc#header x bin)
   flip_vol = 1  [cli]
   voltage = 300.0  [discovered]  (Position_16_3.mdoc#Voltage)
   [ok ] rows value=29 expected=29
   [ok ] dark_sections value=[2, 29]
   WARNING: paths defaulted to 'relative'; set it with --paths or config key series.Position_16_3.paths
   WARNING: dose_convention defaulted to 'exclusive'; set it with --dose-convention or config key series.Position_16_3.dose_convention
wrote cets/at3.cets.json (1 region(s)) + at3.cets-companion.json
report: cets/at3.cets.report.json
```

Output:

```
cets/
├── at3.cets.json              # the CETS dataset: one Region per series (tilt series, movie stacks,
│                              #   alignment "aretomo3", tomograms)
├── at3.cets-companion.json    # what CETS cannot hold: acquisition order, exposure, kV/Cs, dfHand,
│                              #   AlphaOffset/BetaOffset, Thickness, FlipVol, what was dropped
└── at3.cets.report.json       # per series: provenance of every value, gates, warnings, errors
```

More examples:

```bash
# several directories and single files in one dataset
cets-aretomo3 to-cets run_a/ run_b/ extra/TS_99.aln -o cets/all.cets.json --flip-vol 1

# mdocs and tilt stacks live elsewhere
cets-aretomo3 to-cets aretomo3_out/ -o cets/at3.cets.json --mdoc-dir mdoc/ --tilt-stack-dir stacks/

# no stack and no reconstruction next to the .aln: declare what cannot be discovered
cets-aretomo3 to-cets aretomo3_out/ -o cets/at3.cets.json --pix 1.54 --tomo-size 1200

# a local (patch) alignment: keep its global rows only
cets-aretomo3 to-cets aretomo3_out/ -o cets/at3.cets.json --drop-locals
```

## CETS to AreTomo3

```bash
cets-aretomo3 from-cets cets/at3.cets.json -o aretomo3_in/
```

What happens per region:

1. The alignment and the reference tomogram are selected. A region with one of each needs no flags;
   otherwise `--alignment` and `--tomogram` are required.
2. Volume X rotations are checked. An `.aln` cannot hold them, so any non-zero value is refused unless
   `--drop-x-rotation`.
3. The `.aln` is written with NumPatches 0, the centre convention corrected on the way out, and one
   DarkFrame line per tilt image without an alignment.
4. `_TLT.txt` gets the acquisition index and dose from the companion, `--acq-order` or `--dose-per-tilt`.
   Without any of these, a one-column file is written with a warning.
5. `_CTF.txt` is written when every tilt image carries CTF metadata. The defocus hand column is added
   when known (companion or `--defocus-hand`).
6. The written `.aln` and `_TLT.txt` are re-read and checked against the document.

Output:

```
aretomo3_in/
├── TS_01.aln                  # SEC ROT GMAG TX TY SMEAN SFIT SCALE BASE TILT; DarkFrame lines
├── TS_01_TLT.txt              # tilt, acquisition index, dose per raw section
├── TS_01_CTF.txt              # only when every image carries CTF metadata
├── TS_01.mrc -> /data/.../TS_01.mrc    # symlink to the tilt stack when it is local (--no-link-stack)
└── cets_aretomo3.report.json  # provenance, gates and the -Cmd 2 hint per series
```

The report and the console carry the reconstruction command for the written directory. `-AtBin` is
yours to choose; `-CorrCTF 1` is suggested only when a `_CTF.txt` was written:

```
> AreTomo3 -Cmd 2 -Serial 0 -InPrefix aretomo3_in/TS_01.mrc -OutDir aretomo3_in/out/ -Gpu 0 -PixSize 1.54 -kV 300 -Cs 2.7 -AmpContrast 0.07 -VolZ 1200 -AtBin <bin> -FlipVol 1 -Wbp 1 -CorrCTF 1 -SplitSum 0
```

More examples:

```bash
# a document with several alignments per region (e.g. from Warp and from the portal)
cets-aretomo3 from-cets cets/all.cets.json -o aretomo3_in/ --alignment warp --tomogram TS_01_volume

# only some regions, exposure and order supplied by hand
cets-aretomo3 from-cets cets/all.cets.json -o aretomo3_in/ --region TS_01 --region TS_02 \
    --dose-per-tilt 3.0 --acq-order order/TS_01_TLT.txt

# an alignment that carries a volume X rotation (e.g. from Warp's LevelAngleX): drop it explicitly
cets-aretomo3 from-cets cets/warp.cets.json -o aretomo3_in/ --drop-x-rotation
```

## Command reference

### `cets-aretomo3 to-cets`

```
cets-aretomo3 to-cets [OPTIONS] SOURCES...
```

`SOURCES` are `.aln` files or directories (every `*.aln` inside). Duplicates are removed.

| Option | Meaning | Derived from, when absent |
|---|---|---|
| `-o, --output FILE` | CETS dataset JSON to write (required) | |
| `--name TEXT` | dataset name | the output stem |
| `--pix Å` | tilt-image pixel size | stack header, then mdoc `PixelSpacing`; error otherwise |
| `--tomo-size Z` | box thickness in unbinned pixels | `_Vol.mrc` Z extent x bin; error otherwise |
| `--flip-vol 0\|1\|2` | the `-FlipVol` the run used | never inferred; without it `_Vol.mrc` is not described |
| `--mdoc-dir DIR` | where `<stem>.mdoc` lives | next to the `.aln` |
| `--tilt-stack-dir DIR` | where `<stem>.mrc/.st/.mrcs` lives | next to the `.aln` |
| `--dose-convention exclusive\|inclusive` | whether `accumulated_dose` excludes the image's own exposure | `exclusive`, with a warning |
| `--dose-per-tilt D` | constant per-image exposure (e/Å²) | `_TLT.txt`, then mdoc `ExposureDose` |
| `--drop-locals` | keep the global rows of a local alignment | off: local alignments are refused |
| `--no-ctf` | ignore `_CTF.txt` | off |
| `--paths relative\|absolute` | how file paths are written into the document | `relative`, with a warning |
| `--voltage kV`, `--cs mm`, `--amp-contrast F` | companion values only | mdoc voltage; absent otherwise |
| `--fail-fast` | stop at the first failing series | continue, exit 1 at the end |
| `--overwrite` | replace existing outputs | error when outputs exist |
| `--config FILE` | YAML config with overrides (see below) | |

### `cets-aretomo3 from-cets`

```
cets-aretomo3 from-cets [OPTIONS] DOCUMENT
```

`DOCUMENT` is a CETS dataset JSON. Its companion manifest is read when present.

| Option | Meaning | Derived from, when absent |
|---|---|---|
| `-o, --output DIR` | AreTomo3 `-Cmd 2` input directory (required) | |
| `--region ID` | region(s) to convert; repeatable | all |
| `--alignment NAME\|N` | alignment instance name or 0-based index | required when a region has several |
| `--tomogram ID` | reference tomogram id | companion; required when a region has several |
| `--tomo-size Z` | thickness in unbinned pixels | the reference tomogram's extent |
| `--drop-x-rotation` | write an alignment that carries volume X rotations, dropping them | off: refused |
| `--acq-order FILE` | `_TLT.txt`-style file with acquisition index (and dose) per raw section | companion |
| `--dose-per-tilt D` | constant per-image exposure (e/Å²) | companion per-image exposure |
| `--defocus-hand -1\|1` | defocus handedness column of `_CTF.txt` | companion; column omitted otherwise |
| `--link-stack / --no-link-stack` | symlink the tilt stack into the output | on |
| `--no-ctf` | do not write `_CTF.txt` | off |
| `--voltage kV`, `--cs mm`, `--amp-contrast F` | values for the `-Cmd 2` hint | companion; `<?>` in the hint otherwise |
| `--fail-fast`, `--overwrite`, `--config FILE` | as above | |

## Values, defaults and the config file

Every value a command needs resolves through one chain:

```
CLI flag  >  --config FILE  >  companion manifest  >  discovered (headers, mdoc, _TLT.txt)  >  package default
```

Each resolution is printed as `option = value  [source]` and recorded in the report. Falling back to a
package default always prints a warning that names the flag and the config key. Values without a sane
default (pixel size, box thickness) are errors that name the flag.

The config file keeps flag lists short and lets heterogeneous datasets carry per-series values. Keys are
the option names with underscores; unknown keys are errors.

```yaml
cets:                          # every package and command
  voltage: 300
  cs: 2.7
  amp_contrast: 0.07
  paths: relative
cets-aretomo3:
  to-cets: {flip_vol: 1, dose_convention: exclusive, mdoc_dir: mdoc}
  from-cets: {link_stack: true}
series:                        # per tilt series; wins over the command section
  TS_01: {pix: 1.54, tomo_size: 1200}
  TS_02: {pix: 1.34}
```

```bash
cets-aretomo3 to-cets aretomo3_out/ -o cets/at3.cets.json --config cets.yaml
```

## What is refused, what is dropped

- Local (patch) alignments, unless `--drop-locals`. The dropped part is listed in the companion.
- `.aln` files whose SEC / DarkFrame bookkeeping is inconsistent.
- On export, volume X rotations unless `--drop-x-rotation`.
- `_Vol.mrc` written with `-FlipVol 0` (XZY layout) is never described as a tomogram.
- `_CTF.txt` score and resolution columns are not preserved; the writer emits placeholders.

## Development

```bash
pip install -e '.[dev]'
pytest
pre-commit run --all-files     # ruff 0.11.12, ruff-format, mypy 1.8.0
```

The numerical checks of the codec against an independent reference implementation of the AreTomo3
projection model live in cryoet-alignment (`tests/test_cets_goldens.py`).
