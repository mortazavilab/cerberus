You also may want to compute gene triplets based on a the set of isoforms expressed in a given sample. For this analysis, we use [Swan](https://github.com/mortazavilab/swan_vis) to store our sample metadata and expression information:

```python
import swan_vis as swan

sg = swan.SwanGraph()
sg.add_annotation(input.annot)
sg.add_transcriptome(input.gtf, include_isms=True)
sg.add_abundance(input.ab)
sg.add_abundance(input.gene_ab, how='gene')
sg.add_metadata(input.meta)

sg.save_graph(output.prefix)
```

We can use the expression values and metadata stored in the SwanGraph to determine which isoforms are expressed in each sample, and use `cerberus.get_expressed_triplets()` to calculate these values:

```python
import swan_vis as swan
import cerberus

sg = swan.read(swan_file)
ca = cerberus.read('cerberus_annot.h5')

source_name = 'sample_isos'
metadata_column = 'sample'
min_tpm = 1
df = ca.get_expressed_triplets(sg,
                               obs_col=metadata_column,
                               min_tpm=min_tpm,
                               source=source_name)
ca.add_triplets(df)
```
