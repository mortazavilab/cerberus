Replace transcript IDs from your transcriptome in a GTF with cerberus isoform triplets.

```bash
cerberus replace_gtf_ids \
            --h5 {annot} \
            --gtf {gtf} \
            --source {source} \
            --update_ends \
            --collapse \
            -o {update_gtf}
```

**Notes on behavior:**
* Some transcripts from a transcriptome may be assigned the same triplet. Using the `--collapse` flag will deduplicate these in the abundance file and sum up the counts across each transcript with the same isoform triplet.
* Because cerberus uses regions to represent TSSs and TESs rather than single-bp coordinates, the ends of each of your transcripts are not representable in GTF format. We circumvent this by choosing the furthest upstream boundary of the region for the TSS and the furthest downstream boundary for the TES to report in the GTF when the `--update_ends` option is used. If it is not, the original ends from your transcriptome will be used.
