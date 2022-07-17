# cerberus

Cerberus is a set of tools designed to characterize and enhance transcriptome annotations. Currently cerberus can do the following:
* represent transcript start sites (TSSs) and transcript end sites (TESs) as bed regions rather than single base pair ends
* integrate intron chains from multiple transcriptome annotations (GTFs) to create a transcriptome of the union of them all
* integrate TSSs and TESs from multiple GTFs as well as from outside BED sources to create end annotations from the union of them all
* number intron chains, TSSs, and TESs found by their priority in a reference GTF
* use the enhanced intron chain and 5'/3' end sets to annotate an existing GTF transcriptome and to modify the GTF and corresponding abundance matrices to reflect the new naming scheme / identities of the transcripts

# CLI documentation

Cerberus can be run from the command line or used as an API in Python.

## Workflow

### Calling TSS / TES regions from a transcriptome
Create and merge end regions from a transcriptome annotation (GTF) file.

```
Usage: cerberus gtf_to_bed [OPTIONS]

Options:
  --gtf TEXT       GTF file  [required]
  --mode TEXT      Choose tss or tes  [required]
  -o TEXT          Output file name  [required]
  --dist INTEGER   Distance (bp) to extend regions on either side  [default:
                   50]

  --slack INTEGER  Distance allowable for merging regions  [default: 50]
  --help           Show this message and exit.
```

Example calls:
```bash
cerberus gtf_to_bed \
  --mode tss \
  --gtf tests/files/Canx.gtf \
  -o test_output/Canx_tss.bed \
  --dist 50 \
  --slack 50

cerberus gtf_to_bed \
  --mode tes \
  --gtf tests/files/Canx.gtf \
  -o test_output/Canx_tes.bed \
  --dist 50 \
  --slack 50
```

<!-- Calls to generate test files:
```bash
cerberus gtf_to_bed \
  --mode tss \
  --gtf tests/files/Canx.gtf \
  -o tests/files/Canx_tss.bed \
  --dist 50 \
  --slack 50

cerberus gtf_to_bed \
  --mode tes \
  --gtf tests/files/Canx.gtf \
  -o tests/files/Canx_tes.bed \
  --dist 50 \
  --slack 50

cerberus gtf_to_bed \
  --mode tss \
  --gtf tests/files/Canx_1.gtf \
  -o tests/files/Canx_1_tss.bed \
  --dist 50 \
  --slack 50

cerberus gtf_to_bed \
  --mode tss \
  --gtf tests/files/Canx_2.gtf \
  -o tests/files/Canx_2_tss.bed \
  --dist 50 \
  --slack 50
``` -->

### Calling unique intron chains from a transcriptome
Create a tab-separated file detailing unique intron chains present in a
transcriptome annotation (GTF) file.

```
Usage: cerberus gtf_to_ics [OPTIONS]

Options:
  --gtf TEXT  GTF file  [required]
  -o TEXT     Output file name  [required]
  --help      Show this message and exit.
```

Example call:
```bash
cerberus gtf_to_ics \
  --gtf tests/files/Canx.gtf \
  -o test_output/Canx_ics.tsv
```

<!-- Calls to generate test files:
```bash
cerberus gtf_to_ics \
  --gtf tests/files/Canx.gtf \
  -o tests/files/Canx_ics.tsv

cerberus gtf_to_ics \
  --gtf tests/files/Canx_1.gtf \
  -o tests/files/Canx_1_ics.tsv

cerberus gtf_to_ics \
  --gtf tests/files/Canx_2.gtf \
  -o tests/files/Canx_2_ics.tsv
``` -->

### Aggregate end regions from multiple bed files
Create consensus end regions from multiple bed files. The intent is for some
of these files to come from `cerberus gtf_to_bed`.

```
Usage: cerberus agg_ends [OPTIONS]

Options:
  --input TEXT     Path to config file. Each line containsfile path,whether to
                   add ends (True / False),source name  [required]

  --mode TEXT      Choose tss or tes  [required]
  --slack INTEGER  Distance (bp) allowable for merging regions  [default: 20]
  -o TEXT          Output file name  [required]
  --help           Show this message and exit.
```

Example calls:
```bash
cerberus agg_ends \
  --mode tss \
  --input tests/files/Canx_tss_beds.txt \
  -o test_output/Canx_tss_agg.bed

cerberus agg_ends \
  --mode tes \
  --input tests/files/Canx_tes.bed \
  -o test_output/Canx_tes_agg.bed
```

<!-- Calls to generate test files
```bash
cerberus agg_ends \
  --mode tss \
  --input tests/files/Canx_tss_beds.txt \
  -o tests/files/Canx_tss_agg.bed

cerberus agg_ends \
  --mode tes \
  --input tests/files/Canx_tes.bed \
  -o tests/files/Canx_tes_agg.bed
``` -->

### Aggregate intron chains from multiple intron chain files
Create consensus intron chain annotations from multiple intron chain files
(output from `cerberus gtf_to_ics`).

```
Usage: cerberus agg_ics [OPTIONS]

Options:
  --input TEXT  Path to config file. Each line containsfile path,source name
                [required]

  -o TEXT       Output file name  [required]
  --help        Show this message and exit.
```

Example call:
```bash
cerberus agg_ics \
  --input tests/files/Canx_1_ics.tsv,tests/files/Canx_2_ics.tsv \
  -o test_output/Canx_ic_agg.tsv
```

<!-- Calls to generate test files
```bash
cerberus agg_ics \
  --input tests/files/Canx_ics.tsv \
  -o tests/files/Canx_ic_agg.tsv
``` -->

### Create a reference using aggregated intron chains, TSSs, and TESs
Using regions and intron chains from `cerberus agg_ends` and `cerberus agg_ics` respectively, create a `.h5` reference that holds all the information.

```
Usage: cerberus write_reference [OPTIONS]

Options:
  --tss TEXT  TSS bed file output from `agg_ends`  [required]
  --tes TEXT  TES bed file output from `agg_ends`  [required]
  --ics TEXT  IC tsv file output from `agg_ics`  [required]
  -o TEXT     Output .h5 file name
  --help      Show this message and exit.
```

Example call:
```bash
cerberus write_reference \
  --tss test_output/Canx_tss_agg.bed \
  --tes test_output/Canx_tes_agg.bed \
  --ics test_output/Canx_ic_agg.tsv \
  -o test_output/Canx_ref.h5
```

### Annotate an existing GTF transcriptome with the intron chains, TSSs, and TESs from a cerberus reference
Using the reference from `cerberus write_reference`, determine which end regions and intron chain each transcript
in the input GTF uses and output the results to an h5 transcriptome representation.

```
Usage: cerberus annotate_transcriptome [OPTIONS]

Options:
  --gtf TEXT     GTF file  [required]
  --h5 TEXT      cerberus reference from gen_reference  [required]
  --source TEXT  Name of GTF source  [required]
  -o TEXT        Output file name  [required]
  --help         Show this message and exit.
```

Example call:
```bash
cerberus assign-triplets \
  --gtf tests/files/Canx.gtf \
  --h5 test_output/Canx_ref.h5 \
  --source canx \
  -o test_output/Canx_annot.h5
```

<!-- Calls to generate test files:
```bash
cerberus assign-triplets \
  --gtf tests/files/Canx.gtf \
  --ic tests/files/Canx_ic_agg.tsv \
  --tss_bed tests/files/Canx_tss_agg.bed \
  --tes_bed tests/files/Canx_tes_agg.bed \
  -o tests/files/Canx_triplet.h5
``` -->

### Update transcript ids in abundance file
Using the map generated in `cerberus annotate_transcriptome`, update the transcript ids
and transcript names that are used a TALON abundance matrix with the new
triplet versions of the transcript ids / names

```
Usage: cerberus replace_ab_ids [OPTIONS]

Options:
  --h5 TEXT      cerberus reference from gen_reference  [required]
  --ab TEXT      TALON abundance file to replace ids in  [required]
  --source TEXT  name of source in cerberus object to map from  [required]
  --collapse     collapse transcripts with the same triplets
  -o TEXT        Output file name  [required]
  --help         Show this message and exit.
```

Example call:
```bash
cerberus replace_ab_ids \
  --h5 tests/files/Canx_annot.h5 \
  --ab tests/files/Canx_abundance.tsv \
  --source canx \
  --collapse \
  -o test_output/Canx_triplet_updated_abundance.tsv
```

### Update transcript ids in GTF file
Using the map generated in `cerberus annotate_transcriptome`, update the transcript
ids and transcript names that are used in a GTF file with the new triplet versions
of the ids / names

```
Usage: cerberus replace_gtf_ids [OPTIONS]

Options:
  --h5 TEXT      cerberus reference from gen_reference  [required]
  --gtf TEXT     GTF file to replace ids in  [required]
  --source TEXT  name of source in cerberus object to map from  [required]
  --update_ends  Update ends of transcripts with ends from h5
  --collapse     collapse transcripts with the same triplets
  -o TEXT        Output GTF file name  [required]
  --help         Show this message and exit.
```

Example call:
```bash
cerberus replace_gtf_ids \
  --h5 tests/files/Canx_annot.h5
  --gtf tests/files/Canx.gtf \
  --source canx \
  --collapse \
  -o test_output/Canx_triplet_updated.gtf
```

<!-- ## Utilites

### h5 to tsvs
By default as output from `assign-triplets`, cerberus writes a .h5 file with
4 different tables in it corresponding to
* Unique intron chains
* Unique TSS regions in bed format
* Unique TES regions in bed format
* Mapping of transcripts to their corresponding TSS, intron chain, and TES

If you wish to save tsv versions of each of these files for easier viewing,
you can use this utility to convert it.

```
Usage: cerberus h5-to-tsv [OPTIONS]

Options:
  --h5 TEXT     h5 transcriptome file output from cerberus assign-triplets
                [required]
  --opref TEXT  output file prefix  [required]
  --help        Show this message and exit.
```

Example calls:
```bash
cerberus h5-to-tsv \
  --h5 tests/files/Canx_triplet.h5 \
  --opref test_output/Canx
```

<!-- Calls to generate test files:
```bash
cerberus h5-to-tsv \
  --h5 tests/files/Canx_triplet.h5 \
  --opref tests/files/Canx_triplet
``` --> 
