![C/C++ CI](https://github.com/GeoGenetics/unicorn/actions/workflows/c-cpp.yml/badge.svg?branch=unicorn)

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
git clone --recursive https://github.com/GeoGenetics/unicorn.git
cd unicorn
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
  bamstats    Compute per bam statistics.
  taxstats    Compute per taxid statistics.
  reassign    Filter alignments via EM algorithm.
  alnfilt     Filter alignments based on user-defined criteria.
```

## Commands in Detail

### 1. refstats — Per-reference statistics

Compute statistics for each reference sequence.

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

| # | Column | Description |
|---|--------|-------------|
| 1 | Id | Reference name |
| 2 | Length | Reference length |
| 3 | n_alns | Number of alignments |
| 4 | n_reads | Number of reads (≤ n_alns) |
| 5 | m_readl | Mean read length |
| 6 | std_readl | Standard deviation of read length |
| 7 | md_readl | Median read length |
| 8 | mo_readl | Mode read length |
| 9 | readl_min | Minimum read length |
| 10 | readl_max | Maximum read length |
| 11 | m_alnnm | Mean alignment edit distance |
| 12 | m_alnani | Mean alignment ANI |
| 13 | std_alnani | Standard deviation of ANI |
| 14 | md_alnani | Median ANI |
| 15 | n_covbases | Number of covered bases |
| 16 | m_cov | Mean coverage depth |
| 17 | breath_cov | Breadth of coverage |
| 18 | exp_breath | Expected breadth |
| 19 | breath_ratio | Breadth ratio |
| 20 | m_covcovered | Mean coverage of covered positions |
| 21 | std_covcovered | Std dev of coverage of covered positions |
| 22 | evenness_cov | Evenness of coverage |
| 23 | site_density | Site density |
| 24 | entropy | Coverage entropy |
| 25 | gini | Coverage Gini coefficient |
| 26 | n_entropy | Normalised entropy |
| 27 | n_gini | Normalised Gini coefficient |
| 28 | tad80 | Truncated average depth at 80% of coverage mass |
| 29 | mdust | Mean dust score |
| 30 | std_dust | Standard deviation of dust score |

---

### 2. bamstats — Per-BAM statistics

Compute overall statistics for one or more BAM/SAM files.

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

### 3. taxstats — Per-taxid statistics

Compute statistics grouped by taxonomic ID.

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
  --duplicity                  Enable kmer-based duplicity computation.
                               Adds a duplicity column to the output.
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

| # | Column | Description |
|---|--------|-------------|
| 1 | taxid | Taxonomic ID |
| 2 | name | Taxonomic name |
| 3 | num_accessions | Number of reference sequences |
| 4 | total_length | Total reference length |
| 5 | num_alns | Number of alignments |
| 6 | num_reads | Number of reads |
| 7 | mean_readl | Mean read length |
| 8 | stdev_readl | Standard deviation of read length |
| 9 | median_readl | Median read length |
| 10 | mode_readl | Mode read length |
| 11 | readl_min | Minimum read length |
| 12 | readl_max | Maximum read length |
| 13 | mean_alnnm | Mean alignment edit distance |
| 14 | mean_alnani | Mean alignment ANI |
| 15 | stdev_alnani | Standard deviation of ANI |
| 16 | num_covbases | Number of covered bases |
| 17 | mean_cov | Mean coverage depth |
| 18 | breath_cov | Breadth of coverage |
| 19 | exp_breath | Expected breadth |
| 20 | breath_ratio | Breadth ratio |
| 21 | mean_covcovered | Mean coverage of covered positions |
| 22 | site_density | Site density |

With `--duplicity`, a 23rd column is added:

| 23 | duplicity | Fraction of unique k-mers (proxy for genome complexity) |

---

### 4. reassign — EM algorithm filtering

Filter alignments using an Expectation-Maximization algorithm to reassign reads with multiple alignments. Reimplementation of [bamfilter](https://github.com/genomewalker/bam-filter?tab=readme-ov-file#how-the-reassignment-process-works).

> **Requires query-grouped BAM.** Sort with `samtools sort -n input.bam -o query_grouped.bam`.

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
./unicorn reassign -b input.bam --alpha 0.9 --niter 10 -o reassigned.bam
```

---

### 5. alnfilt — Alignment filtering

Filter alignments based on ANI and alignment score criteria.

> **Requires query-grouped BAM.**

```bash
./unicorn alnfilt [options] -b <in.bam>|<in.sam>
Options:
  -b <str>                     Input bam|sam
  -o <str> | --outbam  <str>   Output BAM file [stdout]
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
./unicorn alnfilt -b input.bam --mode ALLTOP --minani 95.0 -o filtered.bam
```

---

## Example workflow

```bash
# 1. Filter ambiguous alignments (requires query-grouped BAM)
./unicorn reassign -b aligned.bam --alpha 0.8 -o reassigned.bam

# 2. Sort by coordinate for taxstats
samtools sort -o sorted.bam reassigned.bam

# 3. Compute per-taxid statistics
./unicorn taxstats \
  -b sorted.bam \
  -a acc2tax.khash \
  -n names.dmp \
  -d nodes.dmp \
  --rank species \
  --minreads 5 \
  --outstat species_stats.txt

# 4. Compute per-reference statistics
./unicorn refstats -b sorted.bam --minreads 5 > reference_stats.txt

# 5. Generate BAM-level summary
./unicorn bamstats -b sorted.bam --outstat bam_summary.txt
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
| bamstats | Any |
| taxstats | Coordinate-sorted (`SO:coordinate`) |
| reassign | Query-grouped (`SO:queryname` or `GO:query`) |
| alnfilt  | Query-grouped (`SO:queryname` or `GO:query`) |

---

## Testing

```bash
make test
```

## Developers

For development information, see [src/README.md](https://github.com/GeoGenetics/unicorn/tree/unicorn/src)

## License

MIT License — see LICENSE file for details.
