|Chromosome|Start|End|Strand|source|Name|
|---|---|---|---|---|---|
|chr1|169794989|169795129|+|v40|ENSG00000000460_1|
|chr1|169795358|169795459|+|encode|ENSG00000000460_2|

* Each BED entry from BED files input to `cerberus agg_ends` will have an entry here.
* Original coordinates are reported in the `Chromosome`, `Start`, `End`, and `Strand` columns.
* `Source` column has source name as defined in `cerberus agg_ends`
* `Name` column includes the final ID for the Cerberus end region that is in the TSS or TES table.
* If the region was not included in the final reference, the `Name` column will be `NaN`.
