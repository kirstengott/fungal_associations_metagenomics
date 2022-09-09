for i in `ls | sed -e "s/\..*$//" | sed -e "s/_.*$//"`; do
    read=`ls $i*sig`
#    echo $read
    for genome in `ls ../../assembly/*sig`
    do
	bn=`basename $genome`
#	echo $bn $i
	 sourmash compare -k 31 --csv ../../../sourmash/genome_raw_read/${i}_${bn}.csv $read $genome
    done
done
