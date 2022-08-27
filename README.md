# fungal_associations_metagenomics


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

for i in `ls data/raw_read/signatures/`; do bn=`echo $i | sed -e "s/\..*$//" | sed -e "s/_.*$//"`; sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-protozoa-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-protozoa-k31.csv; sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-archaea-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-archaea-k31.csv; sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-fungi-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-fungi-k31.csv ; sourmash gather data/raw_read/signatures/$i db/sourmash/genbank-2022.03-viral-k31.zip --threshold-bp 10000 -o sourmash/${bn}_genbank-2022.03-viral-k31.csv; done


## renaming files

rename 's/11463/AttcepgaCombined_FD/' 11463*
rename 's/11465/AptfungaCombined_FD/' 11465*
rename 's/11475/AttcolgardBottom_FD/' 11475*
rename 's/11478/AttcolfgardenTop_FD/' 11478*
rename 's/11489/CypfungaCombined_FD/' 11489*

## aligning reads
bwa mem ../../assembly/AptfungaCombined_FD_2029527003.a.fna AptfungaCombined_FD.1.TCAG.2009_10_03_02_06_45.F3GY48O01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/AptfungaCombined_FD.bam
bwa mem ../../assembly/AttcepgaCombined_FD_2029527004.a.fna  AttcepgaCombined_FD.1.TCAG.2009_10_03_00_06_46.F3GTV4A01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/AttcepgaCombined_FD.bam
bwa mem ../../assembly/AttcolfgardenTop_FD_2029527005.a.fna  AttcolfgardenTop_FD.1.TCAG.2009_10_16_21_47_09.F36MELC01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/AttcolfgardenTop_FD.bam
bwa mem ../../assembly/AttcolgardBottom_FD_2029527006.a.fna AttcolgardBottom_FD.1.TCAG.2009_10_14_07_47_07.F31SS7V01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/AttcolgardBottom_FD.bam
bwa mem ../../assembly/CypfungaCombined_FD_2030936005.a.fna CypfungaCombined_FD.1.TCAG.2009_10_21_06_46_59.F4EQGSL01.fastq.gz | samtools view -@ 5 -b -o ../../aligned_reads/CypfungaCombined_FD.bam



for i in `ls ~/jgi/*/*COG`; do dir=`dirname $i`; bn=`basename $dir`; base=`basename $i`; cp $i ./${bn}_${base}; done

