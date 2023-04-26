Cerberus can aggregate information to support the existence or creation of TSS or TES
regions across multiple bed files. To perform this aggregation, use `cerberus agg_ends`.

```
Usage: cerberus agg_ends [OPTIONS]

Options:
  --input TEXT     Path to config file. Each line contains file path, whether
                   to add ends (True / False), whether to use as a reference,
                   source name  [required]

  --mode TEXT      Choose tss or tes  [required]
  --slack INTEGER  Distance (bp) allowable for merging regions  [default: 20]
  -o TEXT          Output file name  [required]
  --help           Show this message and exit.
```


**Notes on behavior:**
* TSSs / TESs will be numbered from 1...n for each gene based on
  * MANE status, appris_principal status, and basic_set status of TSS / TES (if tags are present in GTF used to call regions)
  * Order of BED files in config file (earlier BED files will yield lower numbers)
* BED files without gene IDs cannot be used to initialize new TSS / TES regions
* After a region is initialized, the boundaries are fixed. Regions that are considered thereafter are simply added as additional forms of support for the already-initialized regions. The reasoning behind this is to avoid "super regions" that grow as new data is added.

**Input BED file format:**

* Generally, BED files follow the [BED file specifications](https://useast.ensembl.org/info/website/upload/bed.html)
* If you have gene ids in your BED file, they must be in the `name` field (4th column), and must be `_` separated with additional text that makes the name unique. Additionally, the `thickStart` (7th) column must  contain the string `cerberus` to indicate to Cerberus how to parse the gene ids. This is how the files output from `cerberus gtf_to_bed` are automatically formatted.

.. _cerberus_bed_format:

.. mdinclude:: cerberus_bed_format.md

**Input config file format:**

* Comma-separated
* No header
* First BED must have `True` add ends setting
* BED files with `True` add ends setting must have gene ids and strand
* Recommended to use `True` add ends setting for GTFs that you are planning to annotate

|BED file path|Reference|Add ends|Source name|
|---|---|---|---|
|path/to/bed/v40_tss.bed|True|True|v40|
|path/to/bed/encode_tss.bed|False|True|encode|


**Output file format:**

.. _cerberus_agg_bed_format:

.. mdinclude:: cerberus_agg_bed_format.md
