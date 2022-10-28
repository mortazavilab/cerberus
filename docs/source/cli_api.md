### cerberus agg_ends
```
Usage: cerberus agg_ends [OPTIONS]

Options:
  --input TEXT     Path to config file. Each line contains file path, whether to
                   add ends (True / False), whether to use as a reference, source name  [required]

  --mode TEXT      Choose tss or tes  [required]
  --slack INTEGER  Distance (bp) allowable for merging regions  [default: 20]
  -o TEXT          Output file name  [required]
  --help           Show this message and exit.
```


### cerberus agg_ics
```
Usage: cerberus agg_ics [OPTIONS]

Options:
  --input TEXT  Path to config file. Each line contains file path,
                whether to use as reference, source name [required]

  -o TEXT       Output file name  [required]
  --help        Show this message and exit.
```

### cerberus annotate_transcriptome

```
Usage: cerberus annotate_transcriptome [OPTIONS]


Options:
  --gtf TEXT          GTF file  [required]
  --h5 TEXT           cerberus reference from gen_reference  [required]
  --source TEXT       Name of GTF source  [required]
  --gene_source TEXT  Source that is already in cerberus object to use gene
                      names from

  -o TEXT             Output file name  [required]
  --help              Show this message and exit.
```

### cerberus gen_reference

```
Usage: cerberus gen_reference [OPTIONS]

Options:
  --ref_gtf TEXT           Path to config file for each GTF, in priority
                           order. Each line contains file path, whether to add
                           ends, whether to use as reference, source name  [required]

  -o TEXT                  Output h5 file name  [required]

  --ref_tss TEXT           Path to config file for each TSS bed file, in
                           priority order. Each line contains file path,
                           whether to add ends, source. Note: these beds will
                           be ordered after the beds from the GTF file.

  --ref_tes TEXT           Path to config file for each TES bed file, in
                           priority order. Each line contains file path,
                           whether to add ends, source. Note: these beds will
                           be ordered after the beds from the GTF file.

  --gtf_tss_dist INTEGER   Distance (bp) to extend TSS regions on either side
                           of each end from the GTFs  [default: 50]

  --gtf_tss_slack INTEGER  Distance allowable for merging TSS regions from
                           each GTF  [default: 50]

  --gtf_tes_dist INTEGER   Distance (bp) to extend TES regions on either side
                           of each end from the GTFs  [default: 50]

  --gtf_tes_slack INTEGER  Distance allowable for merging TSS regions from
                           each GTF  [default: 50]

  --tss_slack INTEGER      Distance (bp) allowable for merging TSS regions
                           [default: 20]

  --tes_slack INTEGER      Distance (bp) allowable for merging TES regions
  --verbosity INTEGER      Verbosity. Higher numbers mean more output.
                           [default: 1]

  --tmp_dir TEXT           Prefix / file path to save temporary files.
                           [default: temp]

  --keep_tmp               Keep intermediate bed and ic files instead of
                           deleting them

  --help                   Show this message and exit.
```

### cerberus gtf_to_bed

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

### cerberus gtf_to_ics

```
Usage: cerberus gtf_to_ics [OPTIONS]

Options:
  --gtf TEXT  GTF file  [required]
  -o TEXT     Output file name  [required]
  --help      Show this message and exit.
```

### cerberus replace_ab_ids

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

### cerberus replace_gtf_ids

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

### cerberus write_reference

```
Usage: cerberus write_reference [OPTIONS]

Options:
  --tss TEXT  TSS bed file output from `agg_ends`  [required]
  --tes TEXT  TES bed file output from `agg_ends`  [required]
  --ics TEXT  IC tsv file output from `agg_ics`  [required]
  -o TEXT     Output .h5 file name
  --help      Show this message and exit.
```
