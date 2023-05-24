In cases where you already have a finalized transcriptome and don't wish to use the Cerberus tools to generate your transcriptome, you have the option to generate summary gene triplets directly from a GTF using a command line tool.

### cerberus gtf_to_triplets

```
Usage: cerberus gtf_to_triplets [OPTIONS]

Options:
  --gtf TEXT            GTF file to replace ids in  [required]
  --source TEXT         name of source in cerberus object to map from
                        [required]

  --gene_id_col TEXT    Attribute name in GTF w/ gene id
  --gene_name_col TEXT  Attribute name in GTF w/ gene name
  -o TEXT               Output file name  [required]
  --help                Show this message and exit.
```

The output file will be a `CerberusAnnotation` h5 object that is directly compatible with Cereberus' downstream plotting functions.
