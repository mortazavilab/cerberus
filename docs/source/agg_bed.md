cerberus can aggregate information to support the existence or creation of TSS or TES
regions across multiple bed files. To perform this aggregation, use `cerberus agg_ends`.

```bash
cerberus agg_ends \
    --input {input.cfg} \
    --mode {mode} \
    --slack {params.slack} \
    -o {output.bed}
```

**Notes on behavior:**
* TSSs / TESs will be numbered from 1...n for each gene based on
  * MANE status, appris_principal status, and basic_set status of TSS / TES (if tags are present in GTF used to call regions)
  * Order of BED files in config file (earlier BED files will yield lower numbers)
* BED files without gene IDs cannot be used to initialize new TSS / TES regions
* After a region is initialized, the boundaries are fixed. Regions that are considered thereafter are simply added as additional forms of support for the already-initialized regions. The reasoning behind this is to avoid "super regions" that grow as new data is added.

**Input bed file format:**

**Input config file format:**

**Output file format:**
