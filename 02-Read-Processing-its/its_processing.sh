#Activate Conda Enviroment
conda activate 2023-USA-Nigeria
#Renaming ITS American kids data
#renaming R1
ls *_R1_*.fastq.gz | sort -n > num_order_R1.txt #get list of file names
paste -d"\t" num_order_R1.txt new_name_R1.txt > name_R1_map.txt #combine list of old names and new names
awk -F'/t' 'system("cp " $1 " " $2)' name_R1_map.txt #rename the files
#renaming R2
sed 's/R1/R2/' new_name_R1.txt > new_name_R2.txt
ls *_R2_*.fastq.gz | sort -n > num_order_R2.txt
paste -d"\t" num_order_R2.txt new_name_R2.txt > name_R2_map.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R2_map.txt
###Follow R script for DADA2 pipeline
#unite database
wget --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie"  https://doi.plutof.ut.ee/doi/10.15156/BIO/2938067/sh_general_release_25.07.2023.tgz
tar -zxvf sh_general_release_25.07.2023.tgz
#making database
mkdir -p ./db
#reformating for database
sed 's/k__Fungi/\tLineage=Root;rootrank;Fungi;domain/g' ./sh_general_release_dynamic_25.07.2023.fasta | sed 's/p__//g' | sed 's/c__/phylum;/g' | sed 's/o__/class;/g' | sed 's/f__/order;/g' | sed 's/g__/family;/g' | sed 's/s__/genus;/g' | sed '/^>/s/$/;species/' > ./db/unite.fasta
perl build_rdp_taxonomy.pl ./db/unite.fasta
mkdir -p ./db/library 
mv ./db/unite.fasta ./db/library/unite.fna
mkdir -p ./db/taxonomy
mv names.dmp nodes.dmp ./db/taxonomy
mv seqid2taxid.map ./db
#building database
kraken2-build --build --db ./db --threads 6
#taxa assignment
cd ~/usa_nigeria/its
kraken2 --db ~/usa_nigeria/its/unite_db/db \
	--threads 7 \
	--use-names \
	--output rep_set.kraken.out rep_set.fa \
	--unclassified-out rep_set.unclassified.out --confidence 0.01
#added to 
#get just ASV, level classified to, and TAX id
python3 combining_ASV_ASV.py

#find mssing ASVs(gotta find msising ASVs on your own :D)
awk -F"\t" '{print $2, $3}' ../rep_set.kraken.out | sed 's/taxid //' > partial_taxa.txt
for i in `cat miss_ASV.txt`; do grep -w $i partial_taxa.txt; done > missing_text.txt
awk '{print $2, $3}' missing_text.txt | sort | uniq | sed 's/(//' | sed 's/)//' | sed 's/ /\t/'> missing_taxa.txt
#cat together seqid2taxid and missing taxa.txt
cat missing_taxa.txt seq2taxid.txt > all2id.txt

#get no repeats for ASVs
sed '1d' ASVs_ASC.tsv | awk '{print $2, $3, $4, $5}' | uniq | awk '!a[$1]++' | sed -e 's/|/\_/g' > uniq_ASVs.txt
awk '{print $4, $1, $2}' uniq_ASVs.txt| sed 's/.*SH/SH/' | sed 's/_re[pf]s_//' | sed 's/singleton_//' | sed 's/ /\t/g' > cleaned_up.txt #get only SH number, ASV, and level

#clean up fasta file
grep ">" sh_general_release_dynamic_25.07.2023.fasta | sed 's/.*SH/SH/' | sed 's/|re[fp]s|/\t/' | sed 's/|refs_singleton|/\t/' > ../SH2taxa.txt
#use python to combine SH names
python3 sh2taxa.py #this gives everything assigned at the species level

#get just the taxa for nonspecies level
grep -v SH uniq_ASVs.txt | grep -v unclassified | awk '{print $4, $1}' > nonspecies_ASVs.txt #get nonspecies and ASV 
awk '{print $1}' nonspecies_ASVs.txt | sort | uniq > nonspects.txt
#get the taxa for non species
awk '{print $2}' SH2taxa.txt | sed 's/[kpcofgs]__/\t/g' | sed 's/;//g' > nonspec2spec.txt
for i in `cat ./nonspects.txt`; do grep -w -m 1 $i nonspec2spec.txt; done > nonspecstaxa.txt
#combine nonspecies and species so that we can then assign taxonomy in python
paste nonspects.txt nonspecstaxa.txt | sed 's/\t//' | sed 's/\t/;/2g'> nonspecies2taxalvl.txt

for i in `cat ./nonspects.txt`; do sed "s/${i}/${i}_blah/2" nonspecies2taxalvl.txt; done | grep blah > blah.txt # remove everything after the assigned taxa level

sed 's/_blah.*//' blah.txt | sort > non_specs.txt
#manually remove anything that isnt assigned at the right level (should have the same number of lines as nonspects.txt) and upload under the same name
awk '$1!=$4{print $1,$4}' data #check to see if columns match
#adding unknown to the end of file so there are 7 taxonomic levels matches
python3 add_unknowns2nonspecs.py non_specs.txt > full_nonspectaxa.txt 
sed 's/ /\t/' nonspecies_ASVs.txt > nonspecs_ASV.txt #format for joining
python3 nonspecs2specs.py

#get all unclassified ASVs together
grep unclassified uniq_ASVs.txt | awk '{print $1, $4}' > unclass_ASVs.txt
python3 add_unknowns2nonspecs.py unclass_ASVs.txt
python3 add_unknowns2nonspecs.py unclass_ASVs.txt | sed 's/Unclassified/Eukaryote_unknown/' > unclass_ASV.txt
#get everything formatted for concating together
sed -i '1d' assigned_nonspec.tsv
awk '{print $3, $4}' assigned_nonspec.tsv > assigned_nonspec.txt
cat assigned_nonspec.txt unclass_ASV.txt | sed 's/;/\t/g' | sed 's/ /\t/' >  unclass_nonspecs.txt
awk '{print $3, $5}' SH_taxa.tsv | sed '1d' | sed 's/[kpcofgs]__/\t/g' | sed 's/;//g' | sed 's/ //' > specs_taxa.txt
cat specs_taxa.txt unclass_nonspecs.txt > its_taxonomy.txt

#final files for phyloseq
its_taxonomy.txt
rep_set_its.fa
sequence_table_its.merged.txt






















#Link to apprioate NCBI database
#get ncbi taxnomy
#mkdir ~/ref_db
#mkdir ~/ref_db/ncbi_taxonomy
#cd ~/ref_db/ncbi_taxonomy
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
#tar -zxvf new_taxdump.tar.gz

#get NCBI taxa
awk -F"\t" '{print $3}' ../rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ~/ref_db/ncbi_taxonomy/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
# check if taxids exist in ranked lineage file
cat taxids | sed 's/$/|/' | while read line; do grep -c -m 1 ^$line rankedlineage_clean | sed "s/0/$line not found/"; done > missing_check
grep "not found" missing_check # should come back with nothing if all taxids found
cat taxids | sed 's/$/|/' | while read line; do grep -m 1 ^$line rankedlineage_clean || echo $line "no lineage" ; done > lineages
grep "no lineage" lineages # should return empty
# remove taxids from file and reverse field order
sed 's/|/\t/' lineages | awk -F"\t" '{print $2}' | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge asv ids and taxonomy
awk '{print $2}' ../rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt



