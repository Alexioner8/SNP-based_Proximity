
for x in SRX23809078 SRX23809076 SRX23809077 SRX23809075 SRX23809083 SRX23809080 SRX23809082 SRX23809084
do
  srr_ids=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$x" | grep -oP '(?<=<Id>)[^<]+' | xargs -I {} curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={}&retmode=xml" | grep -oP '(?<=<PRIMARY_ID>)[^<]+' | grep '\bSRR') 
  srr_symptomatic=$(echo "$srr_ids" | tr '\n' ' ')

  # Create a directory for the current SRX ID
  mkdir -p $x

  for i in $srr_symptomatic
  do
    # Use fastq-dump to output the files to the newly created directory
    ~/sratoolkit.2.11.2-ubuntu64/bin/fastq-dump -O $x $i
  done
done



# nohup bash /mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/download_SRA_symptomatic.sh &