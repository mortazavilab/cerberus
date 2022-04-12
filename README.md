# cerberus

Cerberus is a set of tools designed to use bed regions to represent transcript
start sites (TSSs) and transcript end sites (TESs).

## Workflow

### Calling TSS / TES regions from a transcriptome
Create and merge end regions from a transcriptome annotation (GTF) file.

```
Usage: cerberus gtf-to-bed [OPTIONS]

Options:
  --mode TEXT      Choose tss or tes  [required]
  --gtf TEXT       GTF file  [required]
  -o TEXT          Output file name  [required]
  --dist INTEGER   Distance (bp) to extend regions on either side
                   Default: 50
  --slack INTEGER  Distance allowable for merging regions
                   Default: 50
  --help           Show this message and exit.
```

Example calls:
```bash
cerberus gtf-to-bed \
  --mode tss \
  --gtf tests/files/Canx.gtf \
  -o tests/files/Canx_tss.bed \
  --dist 50 \
  --slack 50

cerberus gtf-to-bed \
  --mode tes \
  --gtf tests/files/Canx.gtf \
  -o tests/files/Canx_tes.bed \
  --dist 50 \
  --slack 50
```

### Calling unique intron chains from a transcriptome
Create a tab-separated file detailing unique intron chains present in a
transcriptome annotation (GTF) file.

```
Usage: cerberus gtf-to-ics [OPTIONS]

Options:
  --gtf TEXT  GTF file  [required]
  -o TEXT     Output file name  [required]
  --help      Show this message and exit.
```

Example call:
```bash
cerberus gtf-to-ics \
  --gtf tests/files/Canx.gtf \
  -o tests/files/Canx_ics.tsv
```

### Aggregate end regions from multiple bed files
Create consensus end regions from multiple bed files.

```
Usage: cerberus agg-ends [OPTIONS]

Options:
  --mode TEXT   Choose tss or tes  [required]
  --input TEXT  Path to file w/ path to BED files on each line or comma-
                separated  list of file paths; ordered by priority  [required]
  -o TEXT       Output file name  [required]
  --help        Show this message and exit.
```

Example calls:
```bash
cerberus agg-ends \
  --mode tss \
  --input tests/files/Canx_tss_beds.txt \
  -o Canx_tss_agg.bed

cerberus agg-ends \
  --mode tes \
  --input tests/files/Canx_tes.bed \
  -o Canx_tes_agg.bed
```

### Compute triplet IDs for a transcriptome
Using the regions from `cerberus agg-ends` combined with the GTF, determine which
end region and intron chain each transcript uses, and number then in ascending
order based on which transcripts are part of the basic set.

```
Usage: cerberus assign-triplets [OPTIONS]

Options:
  --gtf TEXT      GTF of isoforms  [required]
  --tss_bed TEXT  Bed file of TSS regions  [required]
  --tes_bed TEXT  Bed file of TES regions  [required]
  --opref TEXT    Output file prefix to save beds / gtf w/ triplets
                  [required]

  --help          Show this message and exit.
```

Example call:
```bash
cerberus assign-triplets \
  --gtf tests/files/Canx.gtf \
  --tss_bed Canx_tss_agg.bed \
  --tes_bed Canx_tes_agg.bed \
  --opref Canx_triplet
```

### Update transcript ids
Using the map generated in `cerberus assign-triplets`, update the transcript ids
and transcript names that are used a TALON abundance matrix and GTF with the new
triplet versions of the transcript ids / names

```
Usage: cerberus replace-ids [OPTIONS]

Options:
  --map TEXT    transcript ID map from assign_triplets  [required]
  --gtf TEXT    GTF of isoforms
  --ab TEXT     TALON abundance file
  --collapse    collapse transcripts with the same triplets
  --opref TEXT  Output file prefix to save updated gtf / ab  [required]
  --help        Show this message and exit.
```

Example call:
```bash
cerberus replace-ids \
  --map Canx_triplet_tid_map.tsv \
  --gtf tests/files/Canx.gtf \
  --ab tests/files/Canx_abundance.tsv \
  --collapse \
  --opref Canx_triplet
```
