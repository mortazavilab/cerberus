The most simple gene triplet is calculated from all of the isoforms present in a given transcriptome source. You can calculate and add these triplets to your CerberusAnnotation using the `cerberus.get_source_triplets()` function.

In Python, import cerberus, read in your CerberusAnnotation (output from `cerberus annotate_transcriptome`), compute source triplets, and add them to your CerberusAnnotation:

```python
import cerberus

ca = cerberus.read('cerberus_annot.h5')
df = ca.get_source_triplets()
ca.add_triplets(df)
```

After adding these source triplets, you can save your updated CerberusAnnotation in h5 format. The triplets will be saved in the `ca.triplets` DataFrame slot.
```python
ca.write('cerberus_annot_triplets.h5')
```
