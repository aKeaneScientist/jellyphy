#!/bin/bash
start=`date +%s`
jfishstart=`date +%s`
echo running loop to use jellyfish to count kmers in fortyfive fasta files then pass kmer counts to python program to get distance matrix.

arraystrain=( 78 3662 2888 2890 2578 2701 768 568 2433 2791 36 416 31 388 523 18 2489 524 392 2560 2729 2703 2991 1417 2483 2702 3853 2878 2513 2875 2508 543 1495 2644 2898 4020 3000 2559 2693 573 2739 )

mkdir jfishoutput2
mkdir jfishmercounts2

for index in $(seq 0 40);
do
	jellyfish count -m 14 -s 100M -t1 -C -o /home/keanea/algorithm2020/looptest/jfishoutput2/NCYC${arraystrain[index]}mercounts.jf /home/joShare/Ann-Marie/full_dataset/OG41/NCYC${arraystrain[index]}.fasta
	jellyfish dump -c /home/keanea/algorithm2020/looptest/jfishoutput2/NCYC${arraystrain[index]}mercounts.jf_0 > /home/keanea/algorithm2020/looptest/jfishmercounts2/NCYC${arraystrain[index]}merdump.txt
done
jfishend=`date +%s`
echo Jellyfish k-mer counting execution time was `expr $jfishend - $jfishstart` seconds.

python /home/keanea/algorithm2020/looptest/jellyphyfuncs.py --inputDir jfishmercounts2 --outputFile outmtx41_2620.txt > test41time2620.txt

end=`date +%s`
echo Full program time was `expr $end - $start` seconds.
