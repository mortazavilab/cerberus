The final step of reference generation is to write an h5 representation of the TSSs, ICs, and TESs present in the input data as a series of tables in h5 format.

```bash
cerberus write_reference \
    --tss {input.tss} \
    --tes {input.tes} \
    --ics {input.ics} \
    -o {output.ref}
```

**Input TSS / TES file format:**
* Output from `cerberus agg_ends`.

**Input IC file format:**
* Output from `cerberus agg_ics`.

**Output cerberus reference h5 format:**
