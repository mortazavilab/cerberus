You may want to filter the isoforms used to calculate your gene triplets. For example, maybe you want to only consider isoforms that are expressed above a certain TPM level. To compute gene triplets from each gene from any arbitrary list of isoform triplets, use `cerberus.get_subset_triplets()`:

```python
import cerberus
ca = cerberus.read('cerberus_annot.h5')
source_name = 'highly_expressed'
tids = ['Gene[1,1,1]', 'Gene[1,2,1]']
df = ca.get_subset_triplets(tids, source_name)
ca.add_triplets(df)
```
