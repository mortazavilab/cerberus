
To generate BED files for the TSSs and TESs from a GTF file, use `cerberus gtf_to_bed`.

```
Usage: cerberus gtf_to_bed [OPTIONS]

Options:
  --gtf TEXT       GTF file  [required]
  --mode TEXT      Choose tss or tes  [required]
  -o TEXT          Output file name  [required]

  --dist INTEGER   Distance (bp) to extend regions on either side  [default:
                   50]

  --slack INTEGER  Distance allowable for merging regions  [default: 50]
  --help           Show this message and exit.
```

**Output file format:**

.. _cerberus_bed_format:

.. mdinclude:: cerberus_bed_format.md

To generate tab-separated files for the ICs from a GTF file, use `cerberus gtf_to_ics`.

```
Usage: cerberus gtf_to_ics [OPTIONS]


Options:
  --gtf TEXT  GTF file  [required]
  -o TEXT     Output file name  [required]
  --help      Show this message and exit.
```

**Output file format:**

.. _cerberus_ic_format:

.. mdinclude:: cerberus_ic_format.md
