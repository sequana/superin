# Known issues and fixed bugs

## Open

### `get_INDELs` reads the wrong contig for deleted sequences

- **Where:** `sequana_pipelines/superin/tools.py`, `get_INDELs()` (`genome = FastA(fasta_file).sequences[0]`),
  called from `scatter_plots` and `save_plasmid_insertions_in_fastq` in `superin.rules`.
- **Problem:** for the plasmid analysis, `get_fasta` returns the file given in the samplesheet `reference` column
  (e.g. `vibrio_and_pV827.fa`, SI + chromosomes + plasmid). Only the first sequence is used, which is chromosome 1
  and not the plasmid the reads are mapped on. Deleted sequences (`genome[position_ref:position_ref+length]`) are
  therefore taken from the wrong contig.
- **Impact:** probably none on the current outputs, since only `len(sequence)` is used for deletions. It would give
  wrong results if the deleted sequence content was used (e.g. entropy), and the reference position can exceed the
  length of a short first contig.
- **Not caught by tests:** the test references contain only the plasmid sequence.
- **Possible fix:** select the contig matching the read's `rname` (BAM reference name) instead of `sequences[0]`.

## Fixed

### `dot2svg` fails with "missing output .sequana/rulegraph.svg"

- **Symptom:** the pipeline stops in `dot2svg`; `dot` exits without writing the SVG.
- **Cause:** `samplesheet` was stored in `config.yaml` as given on the command line (relative). The `rulegraph`
  snippet runs `snakemake --rulegraph` from the `rulegraph/` sub-directory, where `samplesheet.csv` is not found.
  `rg.dot` and `rulegraph.dot` were then empty, and `dot` on an empty file produces no output.
- **Fix:** `main.py` stores `os.path.abspath(options.samplesheet)`. `test_full` now uses a relative samplesheet
  path so the regression is covered.
- **Existing working directories:** re-run `sequana_superin`, or set `samplesheet:` to an absolute path in
  `.sequana/config.yaml`.
