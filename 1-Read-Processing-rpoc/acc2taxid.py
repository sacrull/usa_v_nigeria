import pandas as pd
df1 = pd.read_table("accessions",delimiter='\t', header= None)
df1.columns =['acc']
print(df1.head())

df2 = pd.read_table("all_taxids2",delimiter='\t', header= None)
df2.columns =['acc', 'TaxID']
print(df2.head())

combined = pd.merge(df1, df2, on ='acc')
inner_join(df2,df1)

print(combined['acc'].unique())

combined.to_csv('acc_taxids.tsv', sep="\t")