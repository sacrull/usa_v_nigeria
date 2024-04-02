import pandas as pd
df1 = pd.read_table("cleaned_up.txt",delimiter='\t', header= None)
df1.columns =['SH', 'ASV', 'Level']
print(df1.head())

df2 = pd.read_table("../SH2taxa.txt",delimiter='\t', header= None)
df2.columns =['SH', 'Taxa']
print(df2.head())

combined = pd.merge(df1, df2, on ='SH')

print(combined['ASV'].unique())

combined.to_csv('SH_taxa.tsv', sep="\t")