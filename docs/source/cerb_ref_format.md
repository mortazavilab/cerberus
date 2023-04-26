* `CerberusAnnotation()` object saved in h5 format.
* Can be read in in Python using

```python
import cerberus
ca = cerberus.read(<ref name>)
```

TSS / TES regions accessible using:
```python
ca.tss
ca.tes
```

.. mdinclude:: cerb_h5_bed_format.md


Because end matching is less precise than intron chain matching, the Cerberus reference h5 also store the mapping between each input TSS / TES and what region in the reference (ie from `ca.tss` or `ca.tes`), if any, the input region matched too.

These maps are stored in:
```python
ca.tss_map
ca.tes_map
```

.. mdinclude:: cerb_h5_end_map.md


Intron chains accessible using:
```python
ca.ic
```

.. mdinclude:: cerb_h5_ic_format.md
