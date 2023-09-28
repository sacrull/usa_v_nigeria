wget -r --no-parent -A '*.1.gff' https://www.homd.org/ftp/genomes/PROKKA/current/gff/
mv *.gff ~/usa_nigeria/rpoc/rpoc_db/test/
wget -r --no-parent -A '*.1.fna' https://www.homd.org/ftp/genomes/PROKKA/current/fna/
mv *.fna ~/usa_nigeria/rpoc/rpoc_db/test/

grep ">" ALL_genomes.fna | sort | sed 's/|.*//' | uniq | sed 's/>//' > gff_files
grep "gene=rpoC" ALL_genomes.gff | sed 's/\t/tab/g' | sed 's/.*ID=//g' | sed 's/;.*//' > rpoC_ids
gffread -x - -g ALL_genomes.fna ALL_genomes.gff | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' > all_genomes.fasta 
for i in `cat rpoC_ids`; do grep -m 1 -A 1 $i all_genomes.fasta; done > rpoc_genes
sed 's/_.*//' rpoC_ids > seqIDs
for i in `cat seqIDs`; do grep -m 1 $i ALL_genomes.fna | sed 's/ .*//' | sed 's/.*|//'; done > acces
for i in `cat ./acces`; do grep -wm 1 $i ../nucl_gb.accession2taxid; done > taxids
cat acces | while read line; do esummary -db nuccore -id $line | xtract -pattern DocumentSummary -element Caption,TaxId,Title; done > taxids2 2>esummary.err

awk '{print $1, $2}' taxids2 > taxid2
awk '{print $1, $2}' taxids3 > taxid3
awk '{print $1, $2}' taxid4 > taxid4
awk '{print $1, $3}' taxids > taxid1
cat taxid1 taxid2 taxid3 | sort | uniq > acces2taxa

#check which ids were not grabbed
awk '{print $1}' acces2taxa > tax_acc
sed 's/\.[123456789]//' acces > accessions
cat accessions tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
cat missed_ids | while read line; do esummary -db nuccore -id $line | xtract -pattern DocumentSummary -element Caption,TaxId,Title; done > taxids5 2>esummary.err

esummary -db nucleotide -id JAKNGA010000001 | xtract -pattern DocumentSummary -element Caption,TaxId

cat ./acces | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId

cat acces | while read line; do esummary -db nuccore -id $line | xtract -pattern DocumentSummary -element Caption,TaxId,Title; done > taxids2 2>esummary.err











grep "CDS" ../test/ALL_genomes.gff | grep 'gene_name "rpoC"' | awk -F"\t" '{print $1, "\t", $9}' | awk -F";" '{print $1}' | sed 's/transcript_id "//' | sed 's/"//' | grep "\."
awk '{print $1}' full_rpoc.ids | awk -F"_" '{print $2}' > accessions
cat accessions | while read line; do esummary -db nuccore -id $line | xtract -pattern DocumentSummary -element Caption,TaxId,Title; done > taxids 2>esummary.err










cat missed_ids | while read p; do 
    echo $p; 
    esearch -db nuccore -query $p < /dev/null | 
    esummary | 
    xtract -pattern DocumentSummary -def "NA" -element Caption,TaxId,Title
    sleep 3s 
done > taxids6 2>esummary.err2













for i in `cat ./gff_files`; do grep -e $i -e Species -A1 ALL_genomes.gff > $i.gff3; done
grep "gene=rpoC" ALL_genomes.gff | sed 's/\t/tab/g' | sed 's/.*ID=//g' | sed 's/;.*//' > rpoC_ids


for i in `cat gff_files`; do grep "gene=rpoC"  | sed 's/\t/tab/g' | sed 's/.*ID=//g' | sed 's/;.*//'; done > rpoC_ids
for i in `cat gff_files`; do gffread -x - -g $i.fna $i.gff | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' > $i.fasta ; done
for i in `cat rpoC_ids`; do grep -A 1 $i all_genomes.fasta; done > rpoc_genes
grep ">" rpoc_genes | sed 's/.1.*/.1/' > acces


grep "CDS" ./ALL_genomes.gff | grep 'gene_name "rpoC"' | awk -F"\t" '{print $1, "\t", $9}' | awk -F";" '{print $1}' | sed 's/transcript_id "//' | sed 's/"//' | grep "\." > full_rpoc.ids