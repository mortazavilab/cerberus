Cerberus can aggregate information to support the existence or creation of intron chains
across multiple datasets. To perform this aggregation, use `cerberus agg_ics`.

```
Usage: cerberus agg_ics [OPTIONS]

Options:

  --input TEXT  Path to config file. Each line contains file path, whether to
                use as a reference, source name  [required]

  -o TEXT       Output file name  [required]
  --help        Show this message and exit.
```

**Notes on behavior:**
* ICs will be numbered from 1...n for each gene based on
  * MANE status, appris_principal status, and basic_set status of IC (if tags are present in GTF used to call intron chains)
  * Order of IC TSV files in config file (earlier TSV files will yield lower numbers)

**Input IC file format:**
* Output from `cerberus gtf_to_ics`.

.. _cerberus_ic_format:

.. mdinclude:: cerberus_ic_format.md

**Input config file format:**

* Comma-separated
* No header
* BED files with `True` add ends setting must have gene ids and strand
* Recommended to use `True` add ends setting for GTFs that you are planning to annotate

|IC file path|Reference|Source name|
|---|---|---|---|
|path/to/bed/v40_tss.bed|True|True|v40|
|path/to/bed/encode_tss.bed|False|True|encode|

**Output file format:**

.. _cerberus_agg_ic_format:

.. mdinclude:: cerberus_agg_ic_format.md
