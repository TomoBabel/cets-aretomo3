# cets-aretomo3

Converters between AreTomo3 tilt-series alignments and the CETS cryo-ET standard, rigid profile
`cets-rigid/0.1` (see `cryoet-alignment/docs/cets.md`). Global (rigid) alignments only; AreTomo3 patch
locals are refused unless dropped explicitly.

```
cets-aretomo3 to-cets aretomo3/ -o out/run.cets.json --flip-vol 1 --mdoc-dir mdoc
cets-aretomo3 from-cets out/run.cets.json -o aretomo3_in/
```

`to-cets` reads, per `.aln`: the tilt stack header (size, pixel size), `_TLT.txt` / mdoc (stage angles,
acquisition order, exposure, frame names, voltage), `_CTF.txt` (per-section CTF), `_Vol.mrc` (reconstruction
dims; described only when `--flip-vol 1|2` is declared — the header cannot tell them apart). It writes one
`Region` per series (tilt series with every raw section, the alignment, the alignment's bin-1 reference
volume and, when declared, the reconstructed tomogram), the companion manifest and a report.

`from-cets` writes `<stem>.aln` (NumPatches 0, DarkFrame lines for images without alignment),
`<stem>_TLT.txt`, `<stem>_CTF.txt` (when every image carries CTF metadata) and prints the `-Cmd 2` line.
A volume X rotation is refused unless `--drop-x-rotation`.

Every value resolves through CLI flag > `--config cets.yaml` > companion > discovered > package default;
a package default always prints a warning naming the flag and config key.

| option | default / derivation |
|---|---|
| `--pix` | tilt stack header › mdoc `PixelSpacing`; required otherwise |
| `--tomo-size Z` | `_Vol.mrc` header × bin › `.aln` `Thickness`; required otherwise |
| `--flip-vol 0\|1\|2` | never inferred; without it `_Vol.mrc` is not described |
| `--dose-per-tilt` | `_TLT.txt` / mdoc per-image exposure |
| `--dose-convention` | `exclusive` (warned) |
| `--voltage --cs --amp-contrast` | mdoc voltage; companion only |
| `--drop-locals`, `--no-ctf`, `--paths relative\|absolute` | off / off / relative |
| from-cets: `--alignment NAME\|N`, `--tomogram ID` | required when a region has several |
| from-cets: `--acq-order FILE`, `--dose-per-tilt`, `--defocus-hand`, `--link-stack` | companion values; 1-column `_TLT.txt` with a warning otherwise |
