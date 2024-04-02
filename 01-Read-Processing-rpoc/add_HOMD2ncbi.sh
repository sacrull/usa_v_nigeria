wget -r --no-parent -A '*.1.gff' https://www.homd.org/ftp/genomes/PROKKA/current/gff/
mv *.gff ~/usa_nigeria/rpoc/rpoc_db/test/
wget -r --no-parent -A '*.1.fna' https://www.homd.org/ftp/genomes/PROKKA/current/fna/
mv *.fna ~/usa_nigeria/rpoc/rpoc_db/test/

grep ">" ALL_genomes.fna | sort | sed 's/|.*//' | uniq | sed 's/>//' > gff_files
grep "gene=rpoC" ALL_genomes.gff | sed 's/\t/tab/g' | sed 's/.*ID=//g' | sed 's/;.*//' > rpoC_ids
sed 's/_.*//' rpoC_ids > seqIDs


wget https://www.homd.org/ftp/genomes/PROKKA/current/SEQID_info.txt
sed 's/ /\t/g' SEQID_info.txt | sed 's/\t.*https/\thttps/'| sed 's/\t.*GCA_/\tGCA_/' | sed 's/_[^_]*//2g' > seqIDs2GCA #format it so it is seqid and GCA
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
sed 's/ /_/g' assembly_summary_genbank.txt | awk '{print $1,$6}' | sed '1d' | sed '1d' | sed 's/ /\t/g' > gca2taxid #get just gca and taxid colum
sed 's/_.*//' rpoC_ids > rpoc_seqs #format rpoc ref seqs with ID
for i in `cat rpoc_seqs`; do grep $i seqIDs2GCA | awk '{print $2}'; done > GCAs #get the GCA for refrence
for i in `cat GCAs`; do grep $i gca2taxid; done > taxid

#check for missing ids
awk '{print $1}' taxid > tax_acc
awk '{print $2}' seqIDs2GCA | sed '1d'> GCA_rpoc 
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
for i in `cat missed_ids`; do grep $i gca2taxid; done > taxid2

cat taxid taxid2 > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids


wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt
for i in `cat missed_ids`; do grep -m 1 $i assembly_summary_genbank_historical.txt | awk '{print $1, $6}'; done > taxid3 #get the GCA for refrence

cat taxid taxid2 taxid3 > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids

for i in `cat missed_ids`; do grep -m 1 $i assembly_summary_genbank.txt | awk '{print $1, $6}'; done > taxid4 #get the GCA for refrence
cat taxid taxid2 taxid3 taxid4 | sed 's/ /\t/' > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids

#get just uniq GCA
sort all_taxids | uniq > GCA_2_taxid
sort rpoc_seqs | uniq | wc -l #sanity check
for i in `cat rpoc_seqs`; do grep $i seqIDs2GCA; done > seqf2GCA_rpoc

#combine each
python3 GCA2taxid.py
sed '1d' acc_taxids.tsv | awk '{print $4}' > taxids
paste rpoc_seqs taxids | sed 's/\t/|kraken:taxid|/' | sed 's/\t/ /' | sed 's/^/>/' > fixed_headers

wget https://www.homd.org/ftp/genomes/PROKKA/current/ffn/ALL_genomes.ffn

#get sequences
gffread -x - -g ALL_genomes.fna ALL_genomes.gff | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' > all_genomes.fasta 
for i in `cat rpoC_ids`; do grep -m 1 -A 1 $i all_genomes.fasta; done > rpoc_genes
grep -v ">" rpoc_genes > seqs
paste fixed_headers seqs | sed 's/\t/\n/' > rpoc_ref.fa

cat rpoc_ref.fa ./HOMD/rpoc_ref.fa > combined_rpoc_ref.fa

mkdir kraken_rpoc
kraken2-build --download-taxonomy --db kraken_rpoc/
kraken2-build --add-to-library combined_rpoc_ref.fa --db kraken_rpoc/
kraken2-build --build --max-db-size 8000000000 --db kraken_rpoc/