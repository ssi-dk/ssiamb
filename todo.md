
# ssiamb — v1 Project TODO Checklist

Use this as a living, end-to-end checklist from empty repo to a tagged `v0.1` release, Bioconda, BioContainer, and Galaxy ToolShed publication.

Legend:  
- [ ] = not started • [~] = in progress • [x] = done

---

## 0) Pre‑flight & Environment
- [x] Confirm Python **3.12+** is available on dev machine(s).
- [x] Create **private test area** for large binaries (minimap2, bwa‑mem2, samtools, bcftools, bbmap, mosdepth) or rely on conda env later.
- [x] Decide initial **admin reference directory** location for local testing and set `$SSIAMB_REF_DIR`.
- [x] Prepare **two tiny bacterial datasets** (reads + assembly/reference) for manual smoke once mapping/calling paths exist.

---

## 1) Repository Bootstrap (Hatchling + Typer Skeleton)
- [x] Initialize repo `ssi-dk/ssiamb` with `src/` layout.
- [x] Add `pyproject.toml` (Hatchling backend, Python ≥3.12, MIT license).
- [x] Add dependencies: `typer[all]`, `rich`, `pyyaml`, `pandas`, `numpy`, `pysam`, `biopython`.
- [x] Create `src/ssiamb/__init__.py` and `src/ssiamb/version.py` (`__version__ = "0.0.0.dev0"`).
- [x] Create `src/ssiamb/cli.py` with Typer app + subcommands: `self`, `ref`, `summarize` (no-op bodies).
- [x] Add console script: `ssiamb = ssiamb.cli:app`.
- [x] **Acceptance:** `pip install -e .` and `ssiamb --help` show global flags + subcommands.

---

## 2) Core Models & Utilities
- [x] `models.py` dataclasses:
  - [x] `Thresholds(dp_min:int, maf_min:float, dp_cap:int)`
  - [x] `SummaryRow` (fields exactly as in spec.md §4.1)
  - [x] `RunPlan` (resolved paths, mode, mapper/caller, thresholds, outputs)
  - [x] Enums `Mode`, `Mapper`, `Caller`, `DepthTool`, `TSVMode`
- [x] `io_utils.py`:
  - [x] Safe **sample name** validator `^[A-Za-z0-9._-]{1,64}$`
  - [x] Sample name inference from R1/R2/VCF/BAM stem with **error** unless `--sample` provided
  - [x] TSV writer with **overwrite/append/fail** and **atomic write** (temp + rename)
  - [x] Simple file lock (best-effort cross‑platform)
  - [x] md5 helper
- [x] `runner.py` stubs: `run_self`, `run_ref`, `run_summarize` producing placeholder `SummaryRow`
- [x] **Acceptance:** `ssiamb self ... --stdout` prints a valid header + row.

---

## 3) Admin Reference Directory & Species Resolution
- [x] `refdir.py` helpers:
  - [x] Resolve ref dir from `$SSIAMB_REF_DIR` or `--ref-dir`
  - [x] `normalize_species_name("Genus species …") -> "Genus_species"` (strip `subsp.`, `complex`, punctuation)
  - [x] Built‑in alias map (e.g., *Shigella sonnei* → **Escherichia coli**, etc.)
  - [x] `resolve_species_to_fasta(ref_dir, species_or_alias) -> Path` (list available species on error)
  - [x] `parse_accession_from_fasta_header(path) -> str|NA`
- [x] Wire `--species` path in `run_ref` and set `ref_label`.
- [x] **Acceptance:** With a fake ref dir, resolution errors list available species keys.

---

## 4) Bracken Parser & Selection
- [x] `bracken.py`:
  - [x] Parse species‑level rows (taxonomy_lvl == "S")
  - [x] Selection algorithm (min_frac=0.70, min_reads=100000; tie‑break: frac → reads → name)
- [x] Wire ref‑mode precedence: `--reference` > `--species` > `--bracken`
- [x] On **no S‑rows** or thresholds unmet, respect `--on-fail` (default `error`) with helpful suggestion.
- [x] **Acceptance:** Unit tests for parser and selector.

---

## 5) Mapping & Index Handling
- [x] `mapping.py`:
  - [x] `ensure_indexes_self(fasta, mapper)` → build next to FASTA if missing (minimap2 `.mmi`, bwa‑mem2 `.0123`, `.amb`, `.ann`, `.pac`, `.bwt.2bit.64`)
  - [x] `map_fastqs(mapper, fasta, r1, r2, sample, threads) -> sorted.bam` with RG `ID/SM/PL`
- [x] `runner.run_self` calls ensure‑indexes (self‑mode) then maps.
- [x] **Acceptance:** `--dry-run` prints planned mapper + index steps.

---

## 6) Denominator via mosdepth
- [x] `depth.py`:
  - [x] Run mosdepth with `--mapq 30 --thresholds 10 --no-per-base`
  - [x] Parse `.mosdepth.summary.txt`
  - [x] Compute `callable_bases`, `genome_length` (≥500 bp contigs), `breadth_10x`, `mean_depth`
  - [x] Expose included contigs set for filtering numerator
- [x] Wire to `run_self`; compute denominator after mapping (and duplicate policy below).
- [x] **Acceptance:** Unit test parsing with synthetic summary file.

---

## 7) Variant Calling
- [x] `calling.py`:
  - [x] BBTools parity: `pileup.sh ploidy=1 minmapq=20 minbaseq=20` → `callvariants.sh ploidy=1 clearfilters=t minallelefraction=0 minavgmapq=20 minavgbasequality=20`
  - [x] bcftools parity‑like: `mpileup -q20 -Q20 -B -a AD,ADF,ADR,DP` → `call -m --ploidy 1`
- [x] Wire stderr capture and exit code **4** on failure.
- [x] **Acceptance:** `--dry-run` shows caller steps.

---

## 8) VCF Normalization & Atomization
- [x] `vcf_ops.py`:
  - [x] `normalize_and_split(vcf_in, ref) -> vcf_out.bgz` (`bcftools norm -f REF -m -both --atomize` + tabix)
- [x] Wire after calling in `run_*`.
- [x] **Acceptance:** dry-run includes normalization steps.

---

## 9) MAF Extraction & Variant Classification
- [x] `vcf_ops.py`:
  - [x] Detect VARIANT_CLASS ∈ {SNV, INS, DEL}; **canonical SNVs only** for primary metric
  - [x] MAF precedence: AD → DP4 (single‑ALT) → AF (multi‑ALT OK); site skip with WARN if none
  - [x] **Per‑site** MAF = `1 − max(allele_fraction)` (multi‑allelic aggregated)
- [x] **Acceptance:** Unit tests for AD/DP4/AF precedence and multi‑alt behavior.

---

## 10) Grid & Primary Count
- [x] `matrix.py` (implemented as part of `vcf_ops.py`):
  - [x] 100×51 grid with `dp_cap` **clipped to 100**
  - [x] `maf_bin = floor(100*MAF)` (0..50)
  - [x] Insert sites (SNV primary only) and build cumulative counts
  - [x] `count_at(dp_min, maf_min)` returns primary `ambiguous_snv_count`
  - [x] TSV writer for wide matrix (rows=1..100; cols=maf_0..maf_50)
  - [x] Gzip compression for matrix output
- [x] **Acceptance:** Unit tests for bin edges and dp_cap clipping.

---

## 11) Contig Filter Alignment
- [x] Ensure **numerator** counts only from contigs **≥500 bp**.
- [x] Use contig include set from mosdepth summary for consistency.
- [x] **Acceptance:** Unit test ensures short contig sites excluded.

---

## 12) Emitters — VCF & BED
- [x] VCF output (when `--emit-vcf`):
  - [x] Emit only records passing thresholds; set `FILTER=PASS`
  - [x] Add INFO: `AMBIG`, `MAF`, `MAF_BIN`, `DP_CAP`, `VARIANT_CLASS`, `ORIG_FILTER` (unless `--require-pass`)
  - [x] Preserve original INFO/FORMAT; left‑normalized and split; bgzip+tabix; **empty-but-valid** OK
- [x] BED output (when `--emit-bed`):
  - [x] 0‑based half‑open; SNV 1‑bp; INS 1‑bp anchor; DEL spans
  - [x] Columns: `chrom start end name score strand sample variant_class ref alt maf dp maf_bin dp_cap`
  - [x] bgzip+tabix; **empty-but-valid** OK
- [x] CLI flags `--emit-vcf` and `--emit-bed` added to `self` and `ref` commands
- [x] Wire emitter functions into `runner.py` for both `self` and `ref` modes
- [x] **Acceptance:** Emitted files sort/index cleanly with tabix; empty outputs valid.

---

## 13) Matrix & Per‑contig Outputs
- [x] `--emit-matrix` → `SAMPLE.variant_matrix.tsv.gz` (wide cumulative grid)
- [x] `--emit-per-contig` → `SAMPLE.per_contig_summary.tsv` with mean depth
- [x] **Acceptance:** Files produced only when flags set.

---

## 14) Summary TSV, MultiQC TSV, Stdout Mode
- [x] Assemble `SummaryRow` with all fields per spec; format rates/rounding.
- [x] Write `ambiguous_summary.tsv` (flat file): default **overwrite**; support append/fail.
- [x] `--stdout`: print only TSV header+row(s); disallow `--emit-*` simultaneously.
- [x] `--emit-multiqc` → `SAMPLE.multiqc.tsv` (curated fields; no runtime_sec).
- [x] **Acceptance:** Overwrite/append modes work; stdout prints valid TSV.

---

## 15) Reuse Inputs & Duplicates
- [x] Enforce **single‑sample** BAM/VCF in v1.
- [x] Precedence: VCF **wins** if both provided; `summarize` requires both.
- [x] Lenient compatibility thresholds:
  - [x] ≥95% length‑weighted contig overlap
  - [x] ≤1% per‑contig length delta
  - [x] ≤2% total length delta
  - [x] Error with code **3** on failure; message lists mismatches; allow `--force-lenient` if implemented later
- [x] Denominator duplicates:
  - [x] If `--bam` provided: if dup flags present → reuse; else run `samtools markdup` transiently for depth path
- [x] **Acceptance:** Mismatch triggers code 3; dup policy observed.

---

## 16) QC Policies & Mapping Rate
- [x] Compute **mapping rate**: primary reads, MAPQ≥20 / primary total (QC‑fail included).
- [x] Warn‑only thresholds:
  - [x] `breadth_10x < 0.80`
  - [x] `callable_bases < 1e6`
  - [x] ref‑mode: map rate < 0.70
- [x] **Acceptance:** WARNs appear; do not fail run.

---

## 17) Provenance (optional)
- [x] Build per‑sample provenance objects with fields from spec §11.
- [x] When `--emit-provenance`, write consolidated `run_provenance.json` (array of objects).
- [x] **Acceptance:** Only written on flag; keys match schema.

---

## 18) Ref‑mode Orchestration
- [x] Finalize precedence in `run_ref`: `--reference` > `--species` > `--bracken`.
- [x] When `--species` used, set `bracken_* = NA` and skip Bracken.
- [x] On Bracken failure, follow `--on-fail` (default `error`); suggest `--species`/`--reference`.
- [x] **Acceptance:** Dry-run clearly shows resolution path and chosen reference.

---

## 19) Summarize Mode
- [ ] Implement `run_summarize`: require `--vcf` and `--bam`.
- [ ] Reuse denominator pipeline from BAM; count with thresholds/grid; optional outputs.
- [ ] **Acceptance:** Errors if either file missing; stdout mode works.

---

## 20) Dry‑run Planner
- [ ] Implement `--dry-run` across subcommands to print planned external calls, thresholds, mapper/caller, resolved reference label, and output plan; **no file writes**.
- [ ] **Acceptance:** Runs exit 0 with readable plan.

---

## 21) Error Handling & Exit Codes
- [ ] Centralize exception → exit code mapping:
  - [ ] 0 success
  - [ ] 1 CLI/input errors (missing args/files, unsafe sample name, invalid manifest)
  - [ ] 2 ref‑mode selection errors
  - [ ] 3 reuse mismatch too severe
  - [ ] 4 external tool failure
  - [ ] 5 QC failure (only if user opts into `--qc-action fail` later)
- [ ] Include helpful suggestions in messages (available species list, thresholds guidance).
- [ ] **Acceptance:** Manual fault injection yields correct codes/messages.

---

## 22) Packaging Polish & Version
- [ ] Add `--version` flag wired to `__version__`.
- [ ] Fill pyproject metadata: homepage/repository, classifiers.
- [ ] README quickstart: install, minimal `self`/`ref` examples, dry‑run, stdout.
- [ ] **Acceptance:** `ssiamb --version` prints version; README renders.

---

## 23) Galaxy Wrapper (ToolShed)
- [ ] Create `galaxy/ssiamb.xml` with params per spec (Mode, inputs, outputs, Advanced, Logging/QC).
- [ ] Add defaults: threads=4, memory=16 GB, walltime=2 h.
- [ ] Autolink Bioconda/BioContainer metadata.
- [ ] Add **tiny self‑mode functional test** (file set TBD by maintainer) so `planemo test` passes.
- [ ] **Acceptance:** `planemo lint` passes; local test passes with placeholder dataset when ready.

---

## 24) CI/CD — PyPI Trusted Publishing
- [ ] GitHub Actions workflow:
  - [ ] Build sdist/wheel on tag `v*`
  - [ ] Use **PyPI Trusted Publishing** (OIDC); project on PyPI accepts trusted publishing
  - [ ] Upload to PyPI on successful build
- [ ] Protect main; add basic lint/test job (optional per spec).
- [ ] **Acceptance:** Dry-run of workflow shows correct steps; publish on `v0.1` tag.

---

## 25) Bioconda Recipe & BioContainer
- [ ] Create Bioconda recipe with PEP 517 build (`pip install .`), run‑deps for external tools.
- [ ] Ensure it produces a BioContainer automatically.
- [ ] Add maintainer handle: `PovilasMat`.
- [ ] **Acceptance:** Local conda build succeeds; PR opened/merged in bioconda‑recipes.

---

## 26) Admin Reference Helper (optional for CLI users)
- [ ] Implement `ssiamb panel download species.txt --out refs/` (writes in‑package `refs/` by default; fail if not writable).
- [ ] Download **best RefSeq** per policy in spec; write both **minimap2** and **bwa‑mem2** indexes.
- [ ] **Acceptance:** Downloaded files exist; names `Genus_species.fna` + indexes.

---

## 27) Minimal Tests & Smoke
- [ ] Unit tests for pure logic:
  - [ ] Species normalization & aliasing
  - [ ] Bracken selection tie‑breaks
  - | [ ] Grid binning & dp_cap clipping
  - [ ] TSV writer header behavior
  - [ ] VCF MAF precedence logic
- [ ] `scripts/dev_smoke.sh` for local quick checks (`--help`, no-op run to placeholder TSV).
- [ ] **Acceptance:** `pytest -q` passes pure tests; smoke script OK.

---

## 28) Documentation (Minimal)
- [ ] Update README with:
  - [ ] Install from PyPI/Conda
  - [ ] Admin ref dir setup & `$SSIAMB_REF_DIR`
  - [ ] Quick examples for `self`, `ref`, `summarize`
  - [ ] Output overview & flags
  - [ ] Error codes table
- [ ] Add `spec.md` and `prompt_plan.md` to repo root.
- [ ] **Acceptance:** Docs reflect current CLI/behavior.

---

## 29) Release
- [ ] Bump version to **v0.1** in `version.py` and `pyproject.toml`.
- [ ] Tag `v0.1` → triggers GH Actions → publish to PyPI.
- [ ] Open Bioconda PR; verify BioContainer build.
- [ ] Publish Galaxy wrapper to ToolShed under owner **pmat**, category **Sequence Analysis**.
- [ ] Announce availability internally; note `$SSIAMB_REF_DIR` setup for Galaxy.

---

## Definition of Done (DoD)
- [ ] `ssiamb self/ref/summarize` run end‑to‑end with real tools where available.
- [ ] `ambiguous_summary.tsv` populated correctly; optional outputs emitted only on flags.
- [ ] `--dry-run` & `--stdout` behave as specified.
- [ ] Exit codes and error messages match spec.
- [ ] v0.1 published to **PyPI**, Bioconda PR opened, Galaxy ToolShed wrapper available.

---

## Nice‑to‑Have (Post‑v1)
- [ ] Add structured JSON logs (`--json-log`).
- [ ] Parallel by samples (`--jobs N`).
- [ ] Plasmid‑aware summaries and alias file on disk.
- [ ] Optional `--fastp` pre‑processing.
- [ ] Masks (`--mask-bed`) and stricter caller‑matched denominators.
