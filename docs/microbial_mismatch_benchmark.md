# Max-mismatch benchmark — microbial gene-sequence probe design

Benchmark of `probe_designer.py` across `--max-mismatches` for the **gene sequence**
input mode of **microbial** probe design, produced by
[benchmark_microbial_mismatches.py](../benchmark/benchmark_microbial_mismatches.py).

| | |
|---|---|
| Run date | 2026-07-28 |
| Runner | `docker exec` into `probe-designer-worker` (GPU pipeline) |
| Result sets | `benchmark/results/microbial_mismatches_{baseline,warm,k30}/` |
| Figures | `mismatch_figure.svg` in each result directory |

---

## Headline findings

1. **The mismatch budget does not affect runtime.** Wall time is flat at
   **≈ 37.4 s** across `--max-mismatches` 0→4. The alignment stage costs the
   same whatever the error budget, so there is no performance argument for
   keeping the budget low.
2. **The mismatch budget is purely a specificity dial.** Raising it 0→4 grows
   retained alignments 61→127 and reclassifies **0→30 of 61** candidate probes
   as off-target — roughly half the panel lost at mm=4.
3. **k-mer length, not mismatch tolerance, decides the final yield.** At the
   default `--kmer-length 18` the pipeline returns **zero** safe probes at
   *every* mismatch value. At `--kmer-length 30` the same runs return
   **60/58/43/36/31**. This is the dominant effect in the whole sweep and it is
   a k-mer sizing problem, not a mismatch problem — see
   [The k=18 collapse](#the-k18-collapse).
4. **All 20 runs completed and every consistency check passed**, including
   bit-identical repeats.

---

## Setup

**Input** — a pasted gene FASTA (`gene_sequences/<job-id>_pasted_gene.fasta`), sequence
`MGYG000000001_1`: **2,580 bp**, GC **28.9 %**. This gene comes from
`MGYG000000001`, a member of the catalog being searched, so the pipeline's
self-match exclusion is exercised.

**Held fixed across the sweep**

| Parameter | Value |
|---|---|
| `--species` | `gut-microbe` |
| `--probe-length` | 36 bp |
| `--tm-range` | 42–47 °C |
| `--gc-range` | 40–80 % |
| `--max-bulges` | 0 — the budget stays purely substitutional |
| alignment targets | microbiome **+** host transcriptome |

**Reference sequence searched — 12.3 GB unique / 13.0 GB scanned**

| Database | Files | Size | |
|---|---:|---:|---|
| `gut-microbe` catalog (`.fa`/`.fna`) | 4,744 | 11.65 GB | |
| Human transcript chunks, current set | 94 | 0.66 GB | |
| Human transcript chunks, stale duplicate | 94 | 0.66 GB | ⚠ redundant |
| **Unique reference** | **4,838** | **12.31 GB** | |
| **Actually scanned** | **4,932** | **12.97 GB** | |

The 4,932 file count matches the pipeline's own log line (`Using 4932
pre-chunked transcript files`). Only sequence is counted: the `.jf` k-mer count
sidecars and `.npz` caches in the same directories are derived indexes, not
searched data, and add ~180 GB if naively included via `du`.

> **⚠ The human transcriptome is duplicated on disk and scanned twice.**
> `data/gencode_data/human/transcript_chunks/` holds two complete parallel chunk
> sets — `…transcripts.filtered.part1–94.fa` (2025-11-18) and
> `…transcripts.part000–093.fa` (2026-05-01, written by
> [split_human_transcripts.py](../split_human_transcripts.py)). They are the same
> data: **279,905 set-identical transcript IDs**, **643,188,897 bases each**, and
> **0 of 279,905 records differ in sequence** — only chunk boundaries and
> filename numbering differ. The `filtered` label is a misnomer; nothing is
> filtered out.
>
> Because `probe_designer.py` globs every `.fa` in that directory, all 188 files
> are loaded and each probe is aligned against the human transcriptome twice.
> This costs ~5 % of alignment time and inflates the host k-mer counts, but does
> **not** affect any conclusion below: `host_rejected` was 0 in all 20 runs and
> the human check reported 0 matches, and duplicate identical data cannot turn a
> zero into a non-zero. The `gut-microbe` catalog is unaffected — all 4,744
> genomes are distinct (the 16 `MGYG*.1.fa` files are accession version suffixes,
> not duplicates).
>
> **Fix:** delete the stale `*.filtered.part*.fa` set. Figures re-rendered after
> that will report 12.3 GB.

**Swept** — `--max-mismatches` ∈ {0, 1, 2, 3, 4}.

**The three result sets**

| Set | Repeats | k-mer | Notes |
|---|---|---|---|
| `baseline` | 1 | 18 | Cold page cache — timings are I/O artefacts, see [Timing](#timing) |
| `warm` | 2 | 18 | `--warmup` + n=2; the timing reference |
| `k30` | 1 | 30 | Identical to `baseline` except `--kmer-length 30` |

---

## Results

### Probe generation is invariant

Every run, at every mismatch value, produced an identical candidate set:

```
2,545 positions checked
  245 pass the GC filter
  179 rejected by Tm, 5 by the homopolymer filter
 → 61 candidate probes
```

`--max-mismatches` is accepted by `core/generate_probes.py` but never consulted
during generation, so this is the expected result and the benchmark asserts it
as a hard check. It confirms the sweep isolates the alignment/safety stages.

### Alignment sensitivity grows with the budget

Warm run (n=2, k=18), mean ± sd:

| max mismatches | Align s | Alignments kept | Off-target probes | Self-only (pre-k-mer safe) |
|---:|---:|---:|---:|---:|
| 0 | 36.5 ± 0.3 | 61 | 0 | 61 |
| 1 | 35.9 ± 0.1 | 66 | 3 | 58 |
| 2 | 36.0 ± 0.1 | 84 | 18 | 43 |
| 3 | 35.7 ± 0.2 | 97 | 25 | 36 |
| 4 | 36.3 ± 0.0 | 127 | 30 | 31 |

Two things worth noting:

- **All 61 probes align at every budget** (`aligned_probes` = 61 throughout,
  `non_aligned_probes` = 0). The budget never changes *whether* a probe aligns
  — only whether its alignments are judged self-only or off-target. The
  pipeline's "non-aligned probes" path is therefore never exercised by this
  input.
- **The loss is non-linear.** mm=0→1 costs 3 probes; mm=1→2 costs 15. The
  jump at mm=2 is where near-identical homologues elsewhere in the catalog
  come into range.

### The k=18 collapse

This is the most consequential result. Same alignment behaviour, opposite
outcome, driven only by `--kmer-length`:

| max mismatches | Self-aligned safe | k-mer unsafe (k=18) | **Safe (k=18)** | k-mer unsafe (k=30) | **Safe (k=30)** |
|---:|---:|---:|---:|---:|---:|
| 0 | 61 | 61 | **0** | 1 | **60** |
| 1 | 58 | 58 | **0** | 0 | **58** |
| 2 | 43 | 43 | **0** | 0 | **43** |
| 3 | 36 | 36 | **0** | 0 | **36** |
| 4 | 31 | 31 | **0** | 0 | **31** |

At k=18 the safety screen rejects **100 % of probes at every mismatch value**.
The pipeline still reports `PIPELINE COMPLETED SUCCESSFULLY` and writes an
empty `safe_probes.fa` — a silent total failure.

The k=18 log shows 251 unique 18-mers from the 61 probes and **553 other-gene
matches**. A back-of-envelope check explains this: with the 11.65 GB microbiome
catalog (~1.1 × 10¹⁰ bp — the human 0.66 GB is screened separately) searched on
both strands (~2.3 × 10¹⁰ positions) against 4¹⁸ ≈ 6.9 × 10¹⁰
possible 18-mers, a *random* 18-mer is expected to occur ≈ 0.33 times by
chance. Each 36 bp probe contributes ~4 unique 18-mers, so on chance alone
~73 % of probes would pick up at least one off-target hit. Add real homology
and 100 % rejection follows. At k=30 the same calculation gives ≈ 2 × 10⁻⁸
expected chance hits per k-mer — negligible — so the single rejected probe at
mm=0 reflects genuine cross-reactivity rather than noise.

**Interpretation:** 18-mers do not carry enough information content to be
specific against a catalog this size. This is a property of the catalog size,
not of the input gene, so it should generalise to any `gut-microbe` run at
k=18. It is the calculation above rather than a direct measurement — a sweep of
`--kmer-length` (20, 22, 25, 30) would locate the actual threshold and is the
obvious follow-up.

### Timing

**Runtime is independent of the mismatch budget.** Warm means range 37.1–37.8 s
with sd ≤ 0.3 s; the k30 set independently reproduces this (37.2–37.7 s).
Alignment accounts for ~96 % of wall time (≈ 36 s of ≈ 37.4 s).

The `baseline` set looks very different and should **not** be read as a trend:

| max mismatches | 0 | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|---:|
| baseline wall s (cold) | 152.4 | 220.7 | 100.2 | 37.5 | 37.7 |
| warm wall s | 37.8 | 37.2 | 37.4 | 37.1 | 37.6 |

These are cold page-cache artefacts from first-touching 4,932 pre-chunked
transcript files. The give-away is that they are not even monotonic — mm=1
(220.7 s) is *slower* than mm=0 (152.4 s) — and that they converge on the warm
value by mm=3 as the cache fills. **Always run with `--warmup`**; the baseline
set is retained only as evidence of the cold-start cost.

`kmer_query_seconds` is 0.00 s throughout. The k-mer counts are computed inside
the fused GPU alignment scan (`fused k-mer counts precomputed`), so the query
is a lookup, not a scan — the safety screen is effectively free.

### The host screen never fired

`host_rejected` is 0 in all 20 runs despite `--align-host` being enabled.
The secondary human-transcript check also reported `Probes matching human
transcripts: 0` every time. For this low-GC bacterial gene that is a plausible
true negative, but it means **these runs do not provide evidence that the host
cross-reactivity path works**. A positive control — a gene with known human
homology — is needed to exercise it.

---

## Consistency checks

All checks passed in all three sets:

| Check | Result |
|---|---|
| all runs completed | PASS — 5 + 10 + 5 = 20 runs |
| candidate probes independent of mismatches | PASS — 61 at every value |
| kept alignments non-decreasing | PASS — 61 → 66 → 84 → 97 → 127 |
| aligned probes non-decreasing | PASS — 61 throughout |
| off-target probes non-decreasing | PASS — 0 → 3 → 18 → 25 → 30 |
| non-aligned probes non-increasing | PASS — 0 throughout |
| safe probes non-increasing | PASS — 60 → 58 → 43 → 36 → 31 (k30) |
| repeat runs deterministic | PASS (warm, n=2); SKIP for n=1 sets |

The determinism check is the strongest result here: across two repeats of all
five mismatch values, `candidate_probes`, `alignments_kept`, `aligned_probes`,
`off_target_probes` and `safe_probes` were identical. The GPU pipeline is
reproducible.

---

## Recommendations

1. **Raise the default `--kmer-length` for large microbiome catalogs, or scale
   it with catalog size.** k=18 is unusable against `gut-microbe`; k=30 works.
   Sweep k ∈ {20, 22, 25, 30} to find the threshold rather than adopting 30 by
   default.
2. **Make a zero-probe outcome loud.** The pipeline currently exits
   `PIPELINE COMPLETED SUCCESSFULLY` after rejecting every probe. It should warn
   when the safety screen rejects ~all candidates and name the likely cause.
3. **Choose the mismatch budget on biology, not cost.** Runtime is flat, so the
   budget is a pure sensitivity/specificity trade. mm=0–1 retains 95–100 % of
   the panel; mm≥2 costs a quarter or more. mm=1 is a reasonable default.
4. **Always pass `--warmup`** when timing.
5. **Add a host-homology positive control** so `--align-host` is actually
   exercised.
6. **Delete the duplicated human transcript chunks** (`*.filtered.part*.fa`) —
   ~0.66 GB of redundant scanning per run. Consider making
   `split_human_transcripts.py`'s cleanup catch legacy naming, or having
   `probe_designer.py` warn when a chunk directory holds more than one chunk
   set, so this cannot recur silently.

---

## Reproducing

```bash
cd benchmark

# warm reference sweep (n=2)
./benchmark_microbial_mismatches.py --repeats 2 --warmup \
    --results-dir results/microbial_mismatches_warm

# k-mer length 30
./benchmark_microbial_mismatches.py --warmup --kmer-length 30 \
    --results-dir results/microbial_mismatches_k30

# redraw a figure without re-running the pipeline
./benchmark_microbial_mismatches.py --plot-only results/microbial_mismatches_k30
```

Each result directory holds `runs.csv` (per-run metrics), `summary.txt`
(table + checks), `results.json` (machine readable), `mismatch_figure.svg`
(three-panel figure) and `logs/` (full pipeline stdout per run).

The script exits non-zero if a strict consistency check fails, so it doubles as
a regression test.

## Related

- [benchmark_microbial_performance.py](../benchmark/benchmark_microbial_performance.py)
  — probe length, mismatch and input-length (1–200 kbp) scaling; results in
  `benchmark/results/microbial_performance/`.
- [tool_methodology_and_outputs.md](tool_methodology_and_outputs.md) — pipeline
  stages and output formats.
- [output_interpretation.rst](output_interpretation.rst) — reading `safe_probes.fa`
  and the k-mer match report.
