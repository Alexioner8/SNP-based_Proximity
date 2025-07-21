
for x in  SRX23809073 SRX23809072 SRX23809086 SRX23809085 SRX23809074 SRX23809081 SRX23809079
do
  srr_ids=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$x" | grep -oP '(?<=<Id>)[^<]+' | xargs -I {} curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={}&retmode=xml" | grep -oP '(?<=<PRIMARY_ID>)[^<]+' | grep '\bSRR') 
  srr_asymptomatic=$(echo "$srr_ids" | tr '\n' ' ')

  # Create a directory for the current SRX ID
  mkdir -p $x

  for i in $srr_asymptomatic
  do
    # Use fastq-dump to output the files to the newly created directory
    ~/sratoolkit.2.11.2-ubuntu64/bin/fastq-dump -O $x $i
  done
done



# nohup bash /mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/download_SRA_asymptomatic.sh &