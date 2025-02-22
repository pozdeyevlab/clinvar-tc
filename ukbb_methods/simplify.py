import polars as pl 

clin = pl.read_csv('variant_summary_pathogenic_Dec_30_24.txt', separator=',', null_values='NA', infer_schema_length=10000)

genes = pl.read_csv('genes.txt', separator='\t', has_header=False)
print(clin)
print(genes)

genes_list = genes['column_1'].to_list()
for g in genes_list:
    print(g)
    d = clin.filter(pl.col('GeneSymbol')==g)
    chr = (d['Chromosome'].unique().to_list())
    print(chr)
    d.write_csv(f'sep2/{chr[0]}_{g}.txt', separator=',')