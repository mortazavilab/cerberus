|Chromosome|Start|Stop|Name|Score|Strand|Sources|Novelty|
|---|---|---|---|---|---|---|---|
|chr1|169794989|169795129|ENSG00000000460_1|.|+|v40,v29,encode|Known|
|chr1|300|400|ENSG00000000460_2|.|+|encode|Novel|

* The `ThickStart` column is used to indicate the sources, defined in the config file, for each region. Sources are comma-separated if there was more than one form of support for the region.
* The `ThickEnd` column is used to indicate the novelty of the region. Regions not supported by one of the sources used as references will be `Novel`, others will be `Known`.
