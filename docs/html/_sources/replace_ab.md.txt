Replace transcript IDs from your transcriptome in an abundance matrix with cerberus isoform triplets.

```bash
cerberus replace_ab_ids \
            --h5 {annot} \
            --ab {ab} \
            --source {source} \
            --collapse \
            -o {update_ab}
```

**Notes on behavior:**
* Some transcripts from a transcriptome may be assigned the same triplet. Using the `--collapse` flag will deduplicate these in the abundance file and sum up the counts across each transcript with the same isoform triplet.

**Input abundance file format:**
* Currently cerberus only works with [TALON](https://github.com/mortazavilab/TALON/) abundance matrices, but we want to support other formats in the future.
