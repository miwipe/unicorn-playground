![C/C++ CI](https://github.com/miwipe/unicorn-playground/actions/workflows/c-cpp.yml/badge.svg?branch=dev)

# Unicorn

Unicorn computes alignment-based statistics from BAM/SAM files for metagenomic analysis.


## Dependencies

Unicorn depends on:
- [htslib](https://github.com/samtools/htslib) for BAM/SAM file handling
- [klib](https://github.com/attractivechaos/klib) (included as submodule)

## Installation

### Standard installation
Make sure htslib is installed and available:

```bash
git clone --recursive https://github.com/miwipe/unicorn-playground.git
cd unicorn-playground
make
```

### Alternative build configurations

If htslib is in a non-standard location, set HTSSRC:

```bash
export HTSSRC=/path/to/htslib/
make
```

`/path/to/htslib` must contain `lib` and `include` and be searchable by the linker at runtime.

For conda environments:

```bash
conda install -c conda-forge -c bioconda enhjoerning
```

Or if you installed htslib from conda and want to compile unicorn yourself:

```bash
export HTSSRC=$CONDA_PREFIX
make
```

## Usage

Unicorn provides five commands for different types of alignment statistics or filtering:

```bash
$ ./unicorn
unicorn 2.4.0 39e362b
        Mar 27 2026 13:11:17
./unicorn command [options] -b <in.bam>|<in.sam>
Commands:
  refstats    Compute per reference statistics.
  alnfilt     Filter alignments based on user-defined criteria.
  taxstats    Compute per taxid statistics.
  reassign    Filter alignments via EM algorithm.
  bamstats    Compute per bam statistics.
```

## Commands in Detail

### 1. refstats — Per-reference statistics

Compute statistics for each reference sequence in the BAM file.

```bash
./unicorn refstats [options] -b <in.bam>|<in.sam>
Options:
  -b <str>                     Input bam|sam [Required]
  -t <int>, --threads <int>    Number of threads [4]
  --outbam  <str>              Output BAM file with filtered alignments.
  --outstat <str>              Output statistics file
  --[FILTER] <PARAM>           Apply filter "FILTER" with parameter "PARAM"
      Available filters:
       - minrefl  <int>        Minimum reference length to consider [0]
       - minreads <int>        Minimum number of reads to consider [1]
       - minalnas <int>        Minimum alignment score [-Inf]
       - maxdust  <int>        Maximum alignment dust score [100]
  --withtid                    Report taxid of reference sequence.
                               Requires --acc2tax, --names and --nodes.
  --names   <str>              Taxonomy nodeid to name mapping file.
  --nodes   <str>              Taxonomy nodeid to parent nodeid mapping file.
  --acc2tax <str>              Accession to taxid mapping file or .khash file.
  --verbose                    Print libunicorn's messages.
  -h                           Print this help message
```

**Basic usage:**
```bash
./unicorn refstats -b input.bam > refstats.txt
```

**With filtering:**
```bash
./unicorn refstats -b input.bam --minreads 10 --minrefl 1000 --outstat refstats.txt
```

**With filtering and filtered BAM output:**
```bash
./unicorn refstats -b input.bam --minreads 10 --minrefl 1000 --outbam filtered.bam > refstats.txt
```

**Output format (30 columns):**

#### Read statistics

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 1 | Id | Reference sequence name | — |
| 2 | Length | Reference length (bp) | — |
| 3 | n_alns | Total number of alignments | Includes multi-mappers |
| 4 | n_reads | Number of unique reads (≤ n_alns) | One read may align to multiple references |
| 5 | m_readl | Mean read length (bp) | For aDNA typically 30–80 bp |
| 6 | std_readl | Standard deviation of read length | Low std dev expected for aDNA |
| 7 | md_readl | Median read length (bp) | Robust to outliers |
| 8 | mo_readl | Mode read length (bp) | Most frequent read length |
| 9 | readl_min | Minimum read length (bp) | — |
| 10 | readl_max | Maximum read length (bp) | — |

#### Alignment quality

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 11 | m_alnnm | Mean alignment edit distance | **Lower = fewer mismatches.** Used alongside ANI |
| 12 | m_alnani | Mean alignment ANI (%) | **Higher = more specific hit.** For modern data >95%; for aDNA >90% typical. Values <90% suggest cross-mapping |
| 13 | std_alnani | Standard deviation of ANI | Low = consistent mapping quality |
| 14 | md_alnani | Median ANI (%) | Robust to outliers; compare with m_alnani to detect skewed distributions |

#### Coverage

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 15 | n_covbases | Number of covered bases | Raw count of positions with ≥1 read |
| 16 | m_cov | Mean coverage depth (×) | Average depth across the entire reference including uncovered positions |
| 17 | breath_cov | Breadth of coverage (0–1) | **Higher = more likely true positive.** Fraction of reference covered by ≥1 read. Low breadth with high depth signals spurious mapping |
| 18 | exp_breath | Expected breadth (0–1) | Theoretical breadth under random uniform read placement given m_cov |
| 19 | breath_ratio | Observed / expected breadth | **Values near 1.0 = coverage consistent with authentic random placement.** Values >>1 suggest pile-ups; values <<1 suggest clustered mapping |
| 20 | m_covcovered | Mean depth of covered positions (×) | Coverage depth restricted to positions that are actually covered |
| 21 | std_covcovered | Std dev of depth on covered positions | **Lower = more even coverage.** High values indicate pile-ups or gaps |
| 22 | evenness_cov | Evenness of coverage (0–1) | **Higher = more even coverage = less likely contaminant** |
| 23 | site_density | Covered bases per kb of reference | Proportion of reference covered, expressed per kilobase |

#### Coverage distribution metrics

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 24 | entropy | Shannon entropy of coverage depth distribution | **Higher = more uniform depth distribution** |
| 25 | gini | Gini coefficient of coverage depth (0–1) | **Lower = more even coverage.** 0 = perfectly uniform; 1 = all reads piled at one position |
| 26 | n_entropy | Normalised Shannon entropy (0–1) | Entropy scaled to [0,1]; comparable across references of different sizes. **Higher = more uniform** |
| 27 | n_gini | Normalised Gini coefficient (0–1) | Gini normalised for reference length. **Lower = more even** |
| 28 | tad80 | Truncated average depth at 80% (×) | Mean depth after removing the top and bottom 10% of coverage mass. More robust than m_cov for references with pile-ups or gaps |

#### Read complexity

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 29 | mdust | Mean DUST score | **Lower = less low-complexity sequence.** High DUST scores indicate repetitive or low-complexity reads |
| 30 | std_dust | Standard deviation of DUST score | High std dev = mixture of complex and low-complexity reads |

---

### 2. alnfilt — Alignment filtering

Filter alignments based on ANI bounds and a user-defined selection strategy. Useful as a pre-processing step before `taxstats`, particularly for eukaryotic data where the EM-based `reassign` is not appropriate.

> **Requires query-grouped BAM.** Sort with `samtools sort -n input.bam -o query_grouped.bam`.

```bash
./unicorn alnfilt [options] -b <in.bam>|<in.sam>
Options:
  -b <str>                     Input bam|sam
  -o <str> | --outbam  <str>   Output BAM file [stdout]
  -t <int> | --threads <int>   Number of threads [4]
  --mode <str>                 Filter mode [ALLTOP]
                               Available modes:
                                ALLTOP  - Keep all best-scoring alignments
                                RNDTOP  - Randomly select one best alignment
                                PCTTOP  - Keep alignments within --pct of best
                                ALL     - Keep all alignments
  --pct <float>                Percentage threshold for PCTTOP mode [0.90]
  --minani <float>             Minimum ANI [90.0]
  --maxani <float>             Maximum ANI [100.0]
  --strictbounds               Remove query if any alignment is out of ANI bounds
  --verbose                    Print libunicorn's messages.
  -h                           Print this help message
```

**Example:**
```bash
./unicorn alnfilt -b query_grouped.bam --mode ALLTOP --minani 95.0 -o filtered.bam
```

---

### 3. taxstats — Per-taxid statistics

Compute statistics grouped by taxonomic ID across all reference sequences assigned to each taxon.

> **Requires coordinate-sorted BAM.** Sort with `samtools sort -o sorted.bam input.bam` before running.

```bash
./unicorn taxstats [options] -b <in.bam>|<in.sam>
Options:
  -b <str>                     Input bam|sam (must be coordinate-sorted)
  -a <str> | --acc2tax <str>   Accession to taxid mapping file or .khash file.
                               Providing a .khash file is much faster.
  -n <str> | --names <str>     Taxonomy names file.
  -d <str> | --nodes <str>     Taxonomy nodes file.
  -t <int> | --threads <int>   Number of threads [4]
  -k <int>                     K-mer size for duplicity computation [17]
  --outstat <str>              Output statistics file [/dev/stdout]
  --[FILTER] <PARAM>           Apply filter "FILTER" with parameter "PARAM"
      Available filters:
       - minrefl  <int>        Minimum reference length [0]
       - minreads <int>        Minimum number of reads per taxid [1]
       - minmani  <float>      Minimum mean ANI per taxid [0]
       - minalnas <int>        Minimum alignment score [-Inf]
       - maxdust  <int>        Maximum alignment dust score [100]
  --duplicity                  Enable duplicity computation (see below).
                               Disabled by default (expensive on large BAMs).
  --filelist <str>             File containing input file paths, one per line.
  --rank <str>                 Taxonomic rank to summarise by [species]
  --verbose                    Print libunicorn's messages.
  -h                           Print this help message
```

**Example:**
```bash
./unicorn taxstats -b sorted.bam \
  -a acc2tax.khash -n names.dmp -d nodes.dmp \
  --rank species --outstat species_stats.txt
```

**With duplicity:**
```bash
./unicorn taxstats -b sorted.bam \
  -a acc2tax.khash -n names.dmp -d nodes.dmp \
  --duplicity --outstat species_stats.txt
```

**Output format:**

Without `--duplicity` (22 columns):

#### Reference and read statistics

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 1 | taxid | NCBI Taxonomic ID | — |
| 2 | name | Taxonomic name at target rank | — |
| 3 | num_accessions | Number of reference sequences in this taxid | More accessions = better genome representation |
| 4 | total_length | Total reference length (bp) | Sum across all accessions |
| 5 | num_alns | Total number of alignments | — |
| 6 | num_reads | Number of unique reads | — |
| 7 | mean_readl | Mean read length (bp) | For aDNA typically 30–80 bp |
| 8 | stdev_readl | Standard deviation of read length | — |
| 9 | median_readl | Median read length (bp) | — |
| 10 | mode_readl | Mode read length (bp) | — |
| 11 | readl_min | Minimum read length (bp) | — |
| 12 | readl_max | Maximum read length (bp) | — |

#### Alignment quality

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 13 | mean_alnnm | Mean alignment edit distance | **Lower = fewer mismatches** |
| 14 | mean_alnani | Mean alignment ANI (%) | **Higher = more specific hit.** For modern data >95%; for aDNA >90% typical |
| 15 | stdev_alnani | Standard deviation of ANI | Low = consistent mapping quality across all accessions |

#### Coverage

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 16 | num_covbases | Total covered bases across all accessions | — |
| 17 | mean_cov | Mean coverage depth (×) | Averaged across all accessions and their full lengths |
| 18 | breath_cov | Breadth of coverage (0–1) | **Higher = more likely true positive.** Fraction of total reference length covered |
| 19 | exp_breath | Expected breadth under random placement (0–1) | Theoretical value given mean_cov |
| 20 | breath_ratio | Observed / expected breadth | **Values near 1.0 = authentic random placement.** Values <<1 suggest pile-ups or highly uneven mapping |
| 21 | mean_covcovered | Mean depth of covered positions (×) | Depth restricted to covered bases only |
| 22 | site_density | Covered bases per kb of reference | Proportion of reference covered, per kilobase |

With `--duplicity`, a 23rd column is added:

#### Duplicity

| # | Column | Description | Interpretation |
|---|--------|-------------|----------------|
| 23 | duplicity | Mean number of times each unique k-mer (k=17 by default) was observed across all reads assigned to this taxid | **Values near 1.0 indicate mostly unique k-mers** — consistent with reads from single-copy genomic regions. **Higher values indicate repetitive sequence or multi-copy elements**, where the same k-mers appear in many reads. Can be used to flag taxa where read evidence is dominated by repetitive or conserved sequence rather than unique genomic content |

---

### 4. reassign — EM algorithm filtering

Filter alignments using an Expectation-Maximization algorithm to probabilistically reassign reads that map to multiple references. Reimplementation of [bamfilter](https://github.com/genomewalker/bam-filter?tab=readme-ov-file#how-the-reassignment-process-works).

> **Requires query-grouped BAM.** Sort with `samtools sort -n input.bam -o query_grouped.bam`.

> **Note:** The EM reassignment algorithm is designed for prokaryotic metagenomic data. It is **not recommended for eukaryotic metagenomic classification**, where genome size variation, repetitive elements, and shared conserved regions make probabilistic read reassignment unreliable. For eukaryotic data, use `alnfilt` with explicit ANI thresholds instead.

```bash
./unicorn reassign [options] -b <in.bam>|<in.sam>
Options:
  -b <str>                     Input bam|sam
  -o <str> | --outbam  <str>   Output BAM file [stdout]
  -t <int> | --threads <int>   Number of threads [4]
  --alpha <float>              Score retention scaling factor (0.0, 1.0] [0.80]
  --niter <int>                Max number of EM iterations [5]
  --scale-type <str>           Subject weight scaling type [LENGTH]
                               Available: NONE, LENGTH, SQRTLEN
  --verbose                    Print libunicorn's messages.
  -h                           Print this help message
```

**Example:**
```bash
./unicorn reassign -b query_grouped.bam --alpha 0.9 --niter 10 -o reassigned.bam
```

---

### 5. bamstats — Per-BAM statistics

Compute overall statistics for one or more BAM/SAM files. Useful as a quick quality check on the BAM file before or after filtering.

```bash
./unicorn bamstats [options] -b <in.bam>|<in.sam>
Options:
  -b <str>         Input bam|sam
  --outstat <str>  Output statistics file
  --filelist <str> File containing input file paths, one per line.
  --printdists     Print distributions of read lengths and ANI.
                   Creates <inputname>.rlen.dists.txt and <inputname>.ani.dists.txt
```

**Example:**
```bash
./unicorn bamstats -b input.bam > bam_summary.txt
```

---

## Example workflows

BAM files are expected in specific sort orders depending on the command used. Starting from an unsorted or coordinate-sorted BAM:

```bash
# Sort by read name (required for alnfilt and reassign)
samtools sort -n -o readname_sorted.bam input.bam

# Sort by coordinate (required for taxstats)
samtools sort -o coordinate_sorted.bam input.bam
```

### Prokaryotic metagenomics

```bash
# 1. Filter ambiguous alignments with EM (requires query-grouped BAM)
./unicorn reassign -b readname_sorted.bam --alpha 0.8 -o reassigned.bam

# 2. Sort by coordinate for taxstats
samtools sort -o coordinate_sorted.bam reassigned.bam

# 3. Compute per-taxid statistics
./unicorn taxstats \
  -b coordinate_sorted.bam \
  -a acc2tax.khash \
  -n names.dmp \
  -d nodes.dmp \
  --rank species \
  --minreads 5 \
  --outstat species_stats.txt
```

### Eukaryotic metagenomics

```bash
# 1. Filter alignments by ANI (EM not recommended for eukaryotes)
./unicorn alnfilt \
  -b readname_sorted.bam \
  --mode ALLTOP \
  --minani 95.0 \
  -o filtered.bam

# 2. Sort by coordinate for taxstats
samtools sort -o coordinate_sorted.bam filtered.bam

# 3. Compute per-taxid statistics with duplicity
./unicorn taxstats \
  -b coordinate_sorted.bam \
  -a acc2tax.khash \
  -n names.dmp \
  -d nodes.dmp \
  --rank species \
  --minreads 5 \
  --duplicity \
  --outstat species_stats.txt
```

### Quality checking

```bash
# Per-reference statistics
./unicorn refstats -b coordinate_sorted.bam --minreads 5 > reference_stats.txt

# BAM-level summary with distributions
./unicorn bamstats -b coordinate_sorted.bam --printdists --outstat bam_summary.txt
```

---

## File formats

### Taxonomy files
- **acc2tax**: Tab-separated file mapping accession IDs to taxonomy IDs
- **names.dmp**: NCBI taxonomy names file
- **nodes.dmp**: NCBI taxonomy nodes file
- **.khash**: Binary format for faster acc2tax lookups — generate with `--dumpacc2tax`

### Input requirements

| Command | BAM sort order required |
|---------|------------------------|
| refstats | Any |
| alnfilt  | Query-grouped (`SO:queryname` or `GO:query`) |
| taxstats | Coordinate-sorted (`SO:coordinate`) |
| reassign | Query-grouped (`SO:queryname` or `GO:query`) |
| bamstats | Any |

---

## Testing

```bash
make test
```

## Developers

For development information, see [src/README.md](https://github.com/miwipe/unicorn-playground/tree/dev/src)

## License

MIT License — see LICENSE file for details.
