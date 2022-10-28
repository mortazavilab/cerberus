cerberus can aggregate information to support the existence or creation of intron chains
across multiple datasets. To perform this aggregation, use `cerberus agg_ics`.

```bash
cerberus agg_ics \
  --input {input.cfg} \
  -o {output.agg_ics}
```

**Notes on behavior:**
* ICs will be numbered from 1...n for each gene based on
  * MANE status, appris_principal status, and basic_set status of IC (if tags are present in GTF used to call intron chains)
  * Order of IC TSV files in config file (earlier TSV files will yield lower numbers)

**Input IC file format:**
* Output from `cerberus gtf_to_ics`.

**Input config file format:**

**Output file format:**
