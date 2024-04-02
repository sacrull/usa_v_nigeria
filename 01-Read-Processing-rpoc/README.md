# rpoC Oral Bacteria Database Build v1.0
Hello!
Here you can find the code that was originally used to build the rpoC db v1.0 and preformmated files of this database for qiime2, kraken2, mothur, dada2.

## Code to build kraken2 rpoC database v1.0
### Notes: seqtk and kraken needs to be installed
### 1. Download bacterial assembly data from NCBI
```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
```
### 2. Download all assemblies with coding sequences annotated
```sh
cat assembly_summary.txt | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GCA_/\t\/GCA_/' | awk '{print $0,$NF}' | sed 's/$/_cds_from_genomic.fna.gz/' | sed 's/\t//g' | sed 's/ //' > query
wget -i query &
# make sure you got them all
grep "pattern" wget-log -B 14 | grep "2023" | grep "ftp:" | awk '{print $3}' > missing #replace 2023 with the current year for you
wget -i missing &
# remove duplicates
rm *fna.gz.1
# concatenate together
cat *fna.gz > all_genomes.fna.gz
gzip -d all_genomes.fna.gz
```
Note: Problem with connections can cause some sequences to not download, which is why you check for those. They have a different error message than bacterial assemblies that do not have a coding sequence. Bacterial assemblies without a coding sequence is one of the reasons the numbers between query and sequences downloaded do not much up.
### 3. Pulling rpoC sequences
```sh
grep "gene=rpoC" all_genomes.fna | grep -v "gamma" | awk '{print $1}' | sed 's/>//' > rpoc.ids
seqtk subseq all_genomes.fna rpoc.ids > rpoc.fna
rm all_genomes.fna
```
NOTE: not pulling gamma subunits -- these are found in chloroplasts and are not the same as rpoC in bacteria. Will have 16106 rpoC loci (this number will potentially go up as more assemblies are added to NCBI). 

### 4. Generate kraken2 database from rpoC sequences
4a. Get the accessions for each sequence
```sh
grep ">" rpoc.fna > headers
sed -i 's/|/_/' headers          
grep -v ">" rpoc.fna > seqs
awk -F"_" '{print $2}' headers > accessions
```
4b. Get taxid from the accessions
```sh
#download and unzip the accession2taxid file
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gzip -d nucl_gb.accession2taxid.gz
#grab the taxid from the accesion number
for i in `cat ../accessions`; do grep -w $i nucl_gb.accession2taxid; done > taxids #this step is very slow
```