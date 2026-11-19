#!/bin/bash

start=$(date +%s)
outdir=$1
# genome=$2
# maskedgenome=$3
threads=$2
echo Creating $outdir output folder
mkdir $outdir
cp $genome $outdir/genome.fa
cp $maskedgenome $outdir/genome.masked.fa
echo Processing genome files
echo 1
javac ChangeFormatMulti.java
java ChangeFormatMulti $outdir/genome.fa $outdir/genome.masked.fa $outdir

echo 2
echo Running TRFinder
chmod +x TRFinder/trf409.linux64
# seqkit split -i $outdir/genome.masked.txt -O $outdir/split   # 按序列拆
# parallel -j 64 './TRFinder/trf409.linux64 {} 2 7 7 80 10 50 500 -f -d -h -ngs > {}.info' ::: $outdir/split/*.txt
# ls -v $outdir/split/*.info | xargs cat > $outdir/genome.masked.info
./TRFinder/trf409.linux64 $outdir/genome.masked.txt 2 7 7 80 10 50 500 -f -d -h -ngs > $outdir/genome.masked.info
echo Finish TRFinder

echo 3
javac MaskTandemRepeatsMulti.java
java MaskTandemRepeatsMulti $outdir/genome.masked.txt $outdir/genome.masked.info $outdir/genome.masked_all.fasta

echo 4
echo Running k-mer counting
chmod +x dsk/dsk
chmod +x dsk/dsk2ascii
./dsk/dsk -verbose 0 -file $outdir/genome.masked_all.fasta -kmer-size 25 -abundance-min 2 -nb-cores $threads -out $outdir/genome_Kmer_2.h5
./dsk/dsk2ascii -verbose 0 -file $outdir/genome_Kmer_2.h5 -nb-cores $threads -out $outdir/genome_Kmer_2.txt
sort -k 1 $outdir/genome_Kmer_2.txt > $outdir/genome_Kmer_2.sorted.txt

echo 5
#python FrequencyLine2.py $outdir/genome_Kmer_2.sorted.txt $outdir/genome.masked_all.fasta $outdir/FrequencyLine.txt 25 $threads
javac FrequencyLineOptimized.java
java FrequencyLineOptimized $outdir/genome_Kmer_2.sorted.txt $outdir/genome.masked_all.fasta $outdir/FrequencyLine.txt 25
echo FrequencyLine fin.

echo 6
javac CombineLineMulti.java
java CombineLineMulti $outdir/FrequencyLine.txt $outdir/genome.txt $outdir/CombineLine.txt
echo CombineLine fin.

echo 7
javac GetAllSCSegsMulti.java
java GetAllSCSegsMulti $outdir/genome.masked_all.fasta $outdir/CombineLine.txt $outdir/genome.txt $outdir/AllSegsOfSCN_CRMasked.fasta $outdir/AllSegsOfSCN_CRRemoved.fasta 25
echo GetAllSCSegs fin.

echo 8
javac PutativeSDsSeparated.java
java PutativeSDsSeparated $threads $outdir/AllSegsOfSCN_CRRemoved.fasta $outdir/AllSegsOfSCN_CRMasked.fasta $outdir
echo PutativeSDsSeparated fin.

echo 9
echo Running Lastz
chmod +x Lastz32.sh
chmod +x /home/fanghuixuan/lastz-distrib/bin/lastz_32
cd $outdir
../Lastz32.sh $threads
wait
cd ..

# cp blast.sh $outdir/blast.sh
# chmod +x blast.sh
# cd $outdir
# ./blast.sh $threads
# wait
# cd ..
echo begin cat
cat $outdir/AllSegsOfSCN_CRRemoved*.info > $outdir/AllSegsOfSCN.info
echo end cat

echo 10
javac LastzPostFast.java
java LastzPostFast $outdir/AllSegsOfSCN_CRMasked.fasta $outdir/AllSegsOfSCN.info $outdir/SCN_LastzResult.txt

echo 11
javac FilterAlignmentsWith500NonCRFast.java
java FilterAlignmentsWith500NonCRFast $outdir/SCN_LastzResult.txt $outdir/SCN_LastzResult_500NonCR.txt 

echo 12
javac SortPairsFast.java
java SortPairsFast $outdir/SCN_LastzResult_500NonCR.txt $outdir/SCN_LastzResult_500NonCR_sortPair.txt 
sort -n -k 5 -k 6 -k 7 -k 8 -k 9 -k 10 -k 11 -k 12 -k 14 $outdir/SCN_LastzResult_500NonCR_sortPair.txt > $outdir/SCN_LastzResult_500NonCR_sortPair.sorted.txt

echo 13
javac FilterRepetitivePairsFast.java
java FilterRepetitivePairsFast $outdir/SCN_LastzResult_500NonCR_sortPair.sorted.txt $outdir/SCN_LastzResult_500NonCRfiltered.txt 

echo 14
javac RecoverFast.java
java RecoverFast $outdir/SCN_LastzResult_500NonCRfiltered.txt $outdir/SCN_LastzResult_500NonCR_FilterPair.txt 

echo 15
javac NewExtendPec50Optimized.java
java NewExtendPec50Optimized $outdir/genome.masked_all.fasta $outdir/genome.txt $outdir/SCN_LastzResult_500NonCR_FilterPair.txt $outdir/SCN_LastzResult_500NonCR_NewExtendPec50.txt

echo 16
javac changeFormatPairwiseOptimized.java
java -Xms8g -Xmx64g -XX:+UseG1GC changeFormatPairwiseOptimized $outdir/genome.masked_all.fasta $outdir/Genome_size_Indexes.txt $outdir/SCN_LastzResult_500NonCR_NewExtendPec50.txt $outdir
echo End pairwise SDs detection. Begin analyze mosaic SDs.

echo 17
javac CopyNewMT.java
java -Xms8g -Xmx64g -XX:+UseG1GC CopyNewMT $outdir/genome.masked_all.fasta $outdir/BG_SDIndexes.fasta $outdir/SCN_LastzResult_500NonCR_NewExtendPec50.txt $outdir/EndpointsLine.txt
sort -n -k 1 -k 2 -k 3 $outdir/BG_SDIndexes.fasta > $outdir/BG_SDIndexes.sorted.fasta

echo 18
javac MergeBGOptimized.java
java MergeBGOptimized $outdir/BG_SDIndexes.sorted.fasta $outdir/BG_MosaicSDs.fasta

echo 19
javac ElementSDsOptimized.java
java -Xms8g -Xmx64g -XX:+UseG1GC ElementSDsOptimized $outdir/genome.masked_all.fasta $outdir/EndpointsLine.txt $outdir/BG_MosaicSDs.fasta $outdir

echo 20
javac testLength.java
java testLength $outdir/genome.masked_all.fasta $outdir/SCN_LastzResult_500NonCR_NewExtendPec50.txt $outdir/length.fasta

echo 21
javac NewElementSDsPairOptimized.java
java NewElementSDsPairOptimized $outdir/genome.masked_all.fasta $outdir/genome.txt $outdir/ClustersLine.txt $outdir/ElementSDs_100.fasta $outdir/length.fasta $outdir/SCN_LastzResult_500NonCR_NewExtendPec50.txt $outdir

echo 22
javac SDblockOptimized.java
java SDblockOptimized $outdir/ElementSDs_100.fasta $outdir/ElementSDs_pairwiseEqual.fasta $outdir

echo 23
javac ElementSDsMultiOptimized.java
java ElementSDsMultiOptimized $outdir/Genome_size_Indexes.txt $outdir/ElementSDs_100.fasta $outdir/SCN_LastzResult_500NonCR_NewExtendPec50.txt $outdir/ElementSDs_LengthAndMulti.fasta

echo 24
javac MosaicSDsBlockcompose.java
java MosaicSDsBlockcompose $outdir/Genome_size_Indexes.txt $outdir/MosaicSDs_SDblockIndexes.txt $outdir/ElementSDs_LengthAndMulti.fasta $outdir/blocks.fasta $outdir/BG_MosaicSDs.fasta


# mkdir $outdir/tmp
# mv $outdir/ElementSDs*.fasta $outdir/tmp/
# mv $outdir/length.fasta $outdir/tmp/
# mv $outdir/Clusters_100.fasta $outdir/tmp/
# mv $outdir/BadMosaicSD_100.fasta $outdir/tmp/
# mv $outdir/BG*.fasta $outdir/tmp/
# mv $outdir/AllSegsOfSCN*.fasta $outdir/tmp/
# mv $outdir/*.fa $outdir/tmp/
# mv $outdir/*.fasta $outdir/tmp/
# mv $outdir/*.info $outdir/tmp/
# mv $outdir/genome*.txt $outdir/tmp/
# mv $outdir/Genome_size_Indexes.txt $outdir/tmp/
# mv $outdir/*.h5 $outdir/tmp/
# mv $outdir/FrequencyLine.txt $outdir/tmp/
# mv $outdir/CombineLine.txt $outdir/tmp/
# mv $outdir/EndpointsLine.txt $outdir/tmp/
# mv $outdir/SCN_LastzResult_500NonCR_NewExtendPec50.txt $outdir/tmp/
# mv $outdir/ClustersLine.txt $outdir/tmp/
# mv $outdir/split/ $outdir/tmp/ 
# #rm -r $outdir/tmp
echo end SDquest analysis
end=$(date +%s)
echo "Total time: $((end-start)) seconds"