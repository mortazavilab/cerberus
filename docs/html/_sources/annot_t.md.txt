Assign each transcript in a GTF file a Cerberus TSS, IC, and TES. You can annotate multiple transcriptomes using the same Cerberus reference.

**Note:** It is recommended that for each transcriptome you annotate, that you
have extracted TSSs, ICs, and TESs from the transcriptome's GTF and added
them to the Cerberus reference. Otherwise you run the risk of encountering triplet
features that are unannotated in the Cerberus reference.

```bash
cerberus annotate_transcriptome \
            --gtf {gtf} \
            --h5 {cerberus ref} \
            --source {source name} \
            -o {cerberus annot}
```

**Output file format:**
