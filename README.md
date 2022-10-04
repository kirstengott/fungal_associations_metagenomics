# fungal_associations_metagenomics


## running quast on assemblies

for i in `ls */*fna | grep -v unassembled | grep -v contigs | grep -v genes`; do outdir=`dirname $i`;  bn=`basename ${i%.fna}`; quast.py -o ${outdir}_${bn} -t 10 -m 300 $i; done

## make uniprot SBT database --- functions.. 

sourmash sketch protein -p k=21,k=31,k=51,scaled=1000,abund -o uniprot_sprot.sig uniprot_sprot.fasta 

sourmash index -k 31 --protein uniprot_sprot uniprot_sprot.sig




## make DNA sketches from reads
for i in `ls`; do trim-low-abund.py -C 3 -Z 18 -V -M 2e9 $i; done

for infile in forward/*.abundtrim; do bn=$(basename ${infile} .fastq.gz.abundtrim); bn1=`echo $bn | sed "s/\..*$//"`; sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge ${bn} -o ${bn}.sig ${infile} reverse/${bn1}*abundtrim; done


for infile in *fastq.gz; do bn=$(basename ${infile} .fastq.gz); sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge ${bn} -o ${bn}.sig ${infile}; done


## make DNA sketches for assemblies
sourmash sketch dna -p scaled=1000,k=21,k=31,k=51 *.fna --name-from-first
for infile in forward/*.abundtrim; do bn=$(basename ${infile} .fastq.gz.abundtrim); bn1=`echo $bn | sed "s/\..*$//"`; sourmash sketch translate -p k=21,k=31,k=51,scaled=1000,abund --merge ${bn} -o ${bn}.protein_sig ${infile} reverse/${bn1}*abundtrim; done

sourmash sketch dna -p scaled=1000,k=21,k=31,k=51 /researchdrive_files/FG_metatranscriptomes/*fastq.gz.abundtrim --name-from-first

for i in `ls /researchdrive_files/FG_metatranscriptomes/*fastq.gz.abundtrim`; do bn=`basename $i`; sourmash sketch translate -p k=21,k=31,k=51,scaled=1000,abund -o $bn.protein_sig $i; done

## getting proteins from Trinity assembly
TransDecoder.LongOrfs -t FG2.trinity_out.Trinity.fna; TransDecoder.LongOrfs -t FG3.trinity_out.Trinity.fna
ls | parallel -j 5 python3 ../../scripts/run_rgi.py {}

python3 ../../scripts/run_diamond.py FG2.trinity_out.Trinity.fna.transdecoder_dir/longest_orfs.pep 5 ../../db/uniprot_sprot.dmnd; python3 ../../scripts/run_diamond.py FG3.trinity_out.Trinity.fna.transdecoder_dir/longest_orfs.pep 5 ../../db/uniprot_sprot.dmnd;
TransDecoder.Predict -t FG2.trinity_out.Trinity.fna --retain_blastp_hits FG2.trinity_out.Trinity.fna.transdecoder_dir/diamond_blastp.tsv; TransDecoder.Predict -t FG3.trinity_out.Trinity.fna --retain_blastp_hits FG3.trinity_out.Trinity.fna.transdecoder_dir/diamond_blastp.tsv

ls assembly/*Trinity.fna | parallel python3 ../scripts/run_infernal.py {}

for infile in *.fastq.gz; do bn=$(basename ${infile} .fastq.gz); sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund -o signatures/${bn}.sig ${infile}; done

## identifying genbank sequences in the reads

for i in `ls data/raw_read/signatures/`; 
	do bn=`echo $i | sed -e "s/\..*$//" | sed -e "s/_.*$//"`; 
	sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-protozoa-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-protozoa-k31.csv; 
	sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-archaea-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-archaea-k31.csv; 
	sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-fungi-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-fungi-k31.csv ; 
	sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-viral-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-viral-k31.csv;
done

ls data/raw_read/signatures/ | grep -v AttbisABBM3 | grep -v AttbisABBM2 | parallel -j 3 python3 scripts/run_sourmash_gather.py data/raw_read/signatures/{} db/sourmash/genbank-2022.03-bacteria-k31.sbt.zip sourmash


## identifying genbank sequences in assemblies 

for i in `ls data/assembly/*sig`; 
	do
	python3 scripts/run_sourmash_gather.py $i db/sourmash/genbank-2022.03-viral-k31.zip sourmash/taxa_genomes & 
	python3 scripts/run_sourmash_gather.py $i db/sourmash/genbank-2022.03-fungi-k31.zip sourmash/taxa_genomes &
	python3 scripts/run_sourmash_gather.py $i db/sourmash/genbank-2022.03-archaea-k31.zip sourmash/taxa_genomes &	
	python3 scripts/run_sourmash_gather.py $i db/sourmash/genbank-2022.03-protozoa-k31.zip sourmash/taxa_genomes &
	python3 scripts/run_sourmash_gather.py $i db/sourmash/genbank-2022.08-plant-k31.sbt.zip sourmash/taxa_genomes
done

ls data/assembly/*sig | parallel -j 3 python3 scripts/run_sourmash_gather.py {} db/sourmash/genbank-2022.03-bacteria-k31.sbt.zip sourmash/taxa_genomes







## renaming files

rename 's/11463/AttcepgaCombined_FD/' 11463*

rename 's/11465/AptfungaCombined_FD/' 11465*

rename 's/11475/AttcolgardBottom_FD/' 11475*

rename 's/11478/AttcolfgardenTop_FD/' 11478*

rename 's/11489/CypfungaCombined_FD/' 11489*

## aligning reads
bwa mem ../../assembly/AptfungaCombined_FD_2029527003.a.fna AptfungaCombined_FD.1.TCAG.2009_10_03_02_06_45.F3GY48O01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/AptfungaCombined_FD.bam

bwa mem ../../assembly/AttcepgaCombined_FD_2029527004.a.fna  AttcepgaCombined_FD.1.TCAG.2009_10_03_00_06_46.F3GTV4A01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/AttcepgaCombined_FD.bam

bwa mem ../../assembly/AttcolfgardenTop_FD_2029527005.a.fna  AttcolfgardenTop_FD.1.TCAG.2009_10_16_21_47_09.F36MELC01.fastq.gz | samtools view 
-@ 5 -b -o ../../aligned_reads/AttcolfgardenTop_FD.bam

bwa mem ../../assembly/AttcolgardBottom_FD_2029527006.a.fna AttcolgardBottom_FD.1.TCAG.2009_10_14_07_47_07.F31SS7V01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/AttcolgardBottom_FD.bam

bwa mem ../../assembly/CypfungaCombined_FD_2030936005.a.fna CypfungaCombined_FD.1.TCAG.2009_10_21_06_46_59.F4EQGSL01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/CypfungaCombined_FD.bam



for i in `ls ~/jgi/*/*COG`; do dir=`dirname $i`; bn=`basename $dir`; base=`basename $i`; cp $i ./${bn}_${base}; done

## Annotate COGS
for i in `ls`; do python3 ../../scripts/run_diamond.py $i 5 ../../db/COG/cog-20.dmnd ../COG/${i%.faa}.dmnd.out; done
blastdbcmd -entry all -db nr -out nr.fa
update_blastdb.pl --decompress nr

for i in `ls`; do ~/scripts/parseBlast.R $i $i.parsed 30 0.001; done

for i in `ls`; do ../..//scripts/parseBlast.R $i $i.parsed 30 0.001; done

for i in `ls`; do ../..//scripts/parseBlast.R $i $i.parsed 30 0.001; done

for infile in *fastq.gz; do bn=$(basename ${infile} .fastq.gz); sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge ${bn} -o ${bn}.sig ${infile}; done

sourmash index -k 31 genbank-2022.03-bacteria-k31 genbank-2022.03-bacteria-k31.zip

for i in `ls data/raw_read/signatures/`; do bn=`echo $i | sed -e "s/\..*$//" | sed -e "s/_.*$//"`; sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.08-plant-k31.sbt.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.08-plant-k31.csv; done

ls data/raw_read/signatures/ | parallel -j 2 python3 scripts/run_sourmash_gather.py data/raw_read/signatures/{} db/sourmash/genbank-2022.03-bacteria-k31.sbt.zip sourmash


## making plant kaiju database 

kaiju-mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o nr_plant_1_kaiju nr_plant_1.fa

kaiju-mkfmi nr_plant_1_kaiju






../../../scripts/run_sourmash_compare_rawr_rawr.sh
../../scripts/run_sourmash_compare_g_g.sh
for i in `ls *fastq.gz | sed -e s/..*$// | sed -e s/_.*$//`; do bwa mem ../assembly/$i*fna $i*fastq.gz | samtools view -@ 3 -b -o ../aligned_reads/${i}.bam& done





## running kaiju on test plant database
../bin/kaiju/bin/kaiju-multi -t ../db/kaiju/nodes.dmp -f ../db/kaiju/nr_plant_1_kaiju.fmi  -i ../data/raw_read/AttbisABBM1_FD_11254.5.198435.CGATGT.fastq.gz,../data/raw_read/AttbisABBM2_FD_11254.5.198435.TTAGGC.fastq.gz,../data/raw_read/AttbisABBM3_FD_11254.6.198438.TGACCA.fastq.gz,../data/raw_read/AttcapACBM1_FD_11231.4.197613.ATGTCA.fastq.gz,../data/raw_read/AttcapACBM2_FD_11231.4.197613.CCGTCC.fastq.gz,../data/raw_read/AttcapACBM3_FD_11231.5.197616.GTAGAG.fastq.gz,../data/raw_read/AttlaeALBM1_FD_11231.5.197616.GTCCGC.fastq.gz,../data/raw_read/AttlaeALBM2_FD_11231.6.197619.GTGAAA.fastq.gz,../data/raw_read/AttlaeALBM3_FD_11231.6.197619.GTGGCC.fastq.gz,../data/raw_read/AttsexASBM1_FD_11231.7.197622.GTTTCG.fastq.gz,../data/raw_read/AttsexASBM2_FD_11231.7.197622.CGTACG.fastq.gz,../data/raw_read/AttsexASBM3_FD_11231.8.197625.GAGTGG.fastq.gz,../data/raw_read/FG1.fastq.gz,../data/raw_read/FG2.fastq.gz,../data/raw_read/FG3.fastq.gz -o AttbisABBM1.plant.out,AttbisABBM2.plant.out,AttbisABBM3.plant.out,AttcapACBM1.plant.out,AttcapACBM2.plant.out,AttcapACBM3.plant.out,AttlaeALBM1.plant.out,AttlaeALBM2.plant.out,AttlaeALBM3.plant.out,AttsexASBM1.plant.out,AttsexASBM2.plant.out,AttsexASBM3.plant.out,FG1.plant.out,FG2.plant.out,FG3.plant.out
rename s/_FD// *FD*
for i in `cat outgroup.txt`; do  sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge $i -o $i/${i}.sig $i/${i}*fna; done
for i in `cat outgroup.txt`; do  sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge $i -o $i/${i}.fna.sig $i/${i}*fna; done
