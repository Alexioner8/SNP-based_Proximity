
for x in SRX23809072 SRX23809073 SRX23809074 SRX23809073 SRX23809075 SRX23809076 SRX23809077 SRX23809078 SRX23809079 SRX23809080 SRX23809082 SRX23809083 SRX23809085 SRX23809086
do
fasta=$(ls /mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/$x)
mkdir ${x}_bam
for file in $fasta
do
# Trimming
java -jar /mnt/raidexttmp/Alejandro/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 /mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${x}/$file \
/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${x}/${file}.fq ILLUMINACLIP:/mnt/raidexttmp/Alejandro/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

rm /mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${x}/$file

/mnt/raidexttmp/Alejandro/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 10 --readFilesIn /mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${x}/${file}.fq --genomeDir /mnt/raidexttmp/Alejandro/GRCh38 --outFileNamePrefix ${x}_bam/${file%.fastq} --outSAMtype BAM SortedByCoordinate
done
done

#nohup bash /mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/STAR_human.sh &

