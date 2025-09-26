

for i in SRR20751889 SRR20751887 SRR20751885
do
~/sratoolkit.2.11.2-ubuntu64/bin/fastq-dump --split-files --gzip $i	

mv ""$i""_1.fastq.gz "$i"_S2_L006_R1_001.fastq.gz
mv "$i"_2.fastq.gz "$i"_S2_L006_R2_001.fastq.gz
mv "$i"_3.fastq.gz "$i"_S2_L006_I1_001.fastq.gz 
done




# cd /mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR --> nohup bash Download_SRA.sh &