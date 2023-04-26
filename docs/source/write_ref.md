The final step of reference generation is to write an h5 representation of the TSSs, ICs, and TESs present in the input data as a series of tables in h5 format.

```
Usage: cerberus write_reference [OPTIONS]

Options:

  --tss TEXT  TSS bed file output from `agg_ends`  [required]
  --tes TEXT  TES bed file output from `agg_ends`  [required]
  --ics TEXT  IC tsv file output from `agg_ics`  [required]
  -o TEXT     Output .h5 file name
  --help      Show this message and exit.
```

**Input TSS / TES file format:**
* Output from `cerberus agg_ends`.

.. _cerberus_agg_bed_format:

.. mdinclude:: cerberus_agg_bed_format.md

**Input IC file format:**
* Output from `cerberus agg_ics`.

.. _cerberus_agg_ic_format:

.. mdinclude:: cerberus_agg_ic_format.md

**Output cerberus reference h5 format:**

.. mdinclude:: cerb_ref_format.md
