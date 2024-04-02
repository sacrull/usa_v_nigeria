####yada dada2
#get the rpoc sequences
cd ~/usa_nigeria
awk '{print $1}' Nigeria_USA_meta.txt | sed '1d' > sample_names.txt
for i in `cat ~/usa_nigeria/sample_names.txt`; do mv "$i"_* ~/usa_nigeria/rpoc/usa_nigeria_rawfastq;  done

#make RPOC database
#download all bacterial assembly 
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
#download all assemblies with coding sequences annotated
cat assembly_summary.txt | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GCA_/\t\/GCA_/' | awk '{print $0,$NF}' | sed 's/$/_cds_from_genomic.fna.gz/' | sed 's/\t//g' | sed 's/ //' > query
wget -i query &
cat assembly_summary.txt | grep "Complete Genome" | awk -F"\t" '{print $20}' | grep GCA_001274535.1 | sed 's/\/GCA_/\t\/GCA_/' | awk '{print $0,$NF}' | sed 's/$/_cds_from_genomic.fna.gz/' | sed 's/\t//g' | sed 's/ //'
# make sure you got them all
grep "pattern" wget-log.1 -B 14 | grep "2023" | grep "ftp:" | awk '{print $3}' > missing  #replace 2023 with whatever year you are in
wget -i missing &
# remove duplicates
rm *fna.gz.1
ls | grep .fna.fz >fnas
mkdir fna
for i in `cat fnas`; do mv $i fna ; done
# concatenate together
cat *fna.gz > all_genomes.fna.gz
find . -type f -name "*.fna.gz" -print0 | xargs -0 cat > all_genomes.fna.gz
gzip -d all_genomes.fna.gz
#pull rpoc gene
grep "gene=rpoC" all_genomes.fna | grep -v "gamma" | awk '{print $1}' | sed 's/>//' > rpoc.ids
~/bin/seqtk/seqtk/seqtk subseq all_genomes.fna rpoc.ids > rpoc.fna
rm all_genomes.fna
#generate kraken2 db
#format fasta file
grep ">" rpoc.fna > headers
sed -i 's/|/_/' headers          
grep -v ">" rpoc.fna > seqs
awk -F"_" '{print $2}' headers > accessions
#use asscesions to get taxids
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gzip -d nucl_gb.accession2taxid.gz
#get taxids
for i in `cat ../accessions`; do grep -wm 1 $i nucl_gb.accession2taxid; done > taxids
#get missed taxids
awk '{print $2}' taxids > tax_acc
cat accessions tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
#may have to manually create file with accession \t taxid
awk '{print $2, $3}' taxids > accs_taxids
cat accs_taxids accs_missed > all_taxids #combine misisng taxids with other taxids
sed 's/ /\t/' all_taxids > temp
mv temp all_taxids 
#check all taxids are there
awk '{print $1}' all_taxids > tax_acc
cat accessions tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids #this should be empty now
sort all_taxids | uniq > all_taxids2
#fix headers
awk '{print $2, $3}' acc_taxids.tsv | sed '1d'| sed 's/^/>/' | sed 's/ /|kraken:taxid|/' | sed 's/\t/ /' > fixed_headers
paste fixed_headers ./fna/seqs | sed 's/\t/\n/' > rpoc_ref.fa
#generate database
mkdr kraken_rpoc
kraken2-build --download-taxonomy --db kraken_rpoc/
kraken2-build --add-to-library ./rpoc_ref.fa --db kraken_rpoc
kraken2-build --build --max-db-size 8000000000 --db kraken_rpoc

#taxonomic assignment
kraken2 --db ~/usa_nigeria/rpoc/rpoc_db/kraken_rpoc \
	--threads 6 \
	--use-names \
	--output rep_set.kraken.out rep_set_rpoc.fa \
	--unclassified-out rep_set.unclassified.out --confidence 0.01

#get full taxonomy
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ~/ref_db/ncbi_taxonomy/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
# check if taxids exist in ranked lineage file
cat taxids | sed 's/$/|/' | while read line; do grep -c -m 1 ^$line rankedlineage_clean | sed "s/0/$line not found/"; done > missing_check
grep "not found" missing_check # should come back with nothing if all taxids found
cat taxids | sed 's/$/|/' | while read line; do grep -m 1 ^$line rankedlineage_clean || echo $line "no lineage" ; done > lineages #use the || pipe
grep "no lineage" lineages # should return empty

# remove taxids from file and reverse field order
sed 's/|/\t/' lineages | awk -F"\t" '{print $2}' | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge asv ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt

#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set_rpoc.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt

#make krona plot
ktImportTaxonomy -t 5 -m 3 -o krona.html rep_set.kraken.out 