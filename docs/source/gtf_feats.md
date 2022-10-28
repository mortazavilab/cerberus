
To generate BED files for the TSSs and TESs from a GTF file, use `cerberus gtf_to_bed`.

```bash
cerberus gtf_to_bed \
    --gtf {gtf} \
    --mode {wildcards.mode} \
    --dist {params.dist} \
    --slack {params.slack} \
    -o {output.bed}
```

**Output file format:**


To generate tab-separated files for the ICs from a GTF file, use `cerberus gtf_to_ics`.

```bash
cerberus gtf_to_ics \
    --gtf {input.gtf} \
    -o {output.ics}
```

**Output file format:**
