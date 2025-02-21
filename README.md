# CSON_GENOME
 This is a walk through of how I did the CSON Genome
 
 
 ## Adapter Removal
 
 apptainer version: 1.3.1
 pdadaoterfilt version: Version 3
 
 ```
 apptainer build pbadapterfilt.sif docker://dmolik/pbadapterfilt
 apptainer exec ./pbadapterfilt.sif pbadapterfilt.sh cson_f_hifi.bam -t 16
 gunzip -k cson_f_hifi.filt.fastq.gz
 awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' cson_f_hifi.filt.fastq > cson_f_hifi.filt.fasta
 ```
 
## Mitochondrial Genome
 
 MitoHiFi Version: 3.2.2
 
```
apptainer build mitohifi.sif docker://ghcr.io/marcelauliano/mitohifi:master
apptainer exec ./mitohifi.sif findMitoReference.py --species "Culicoides sonorensis" --outfolder ./mitoref --min_length 14000
apptainer exec ./mitohifi.sif mitohifi.py -r cson_f_hifi.filt.fasta -f ./mitoref/BK065013.1.fasta -g ./mitoref/BK065013.1.gb -t 8 -o 5
```

## Removing Mitochondrial Reads

Minimap2 Version: 2.26
Samtools Version: 1.17

```
apptainer build minimap2_2.26.sif docker://ebird013/minimap2:2.26
apptainer build samtools_1.17.sif docker://ebird013/samtools:1.17
apptainer exec ./minimap2_2.26.sif minimap2 -ax map-hifi final_mitogenome.fasta cson_f_hifi.filt.fastq.gz > mito_f.sam
apptainer exec ./samtools_1.17.sif samtools sort mito_f.sam | samtools view -f 4 | samtools fastq -1 cson_f_hifi.filt.mitorm.fastq.gz
```

## Assembling with HiFi-ASM

HiFi-ASM Version: 0.24.0

```
apptainer build hifiasm_0.24.0.sif docker://quay.io/biocontainers/hifiasm:0.24.0--h5ca1c30_0
apptainer exec ./hifiasm_0.24.0.sif hifiasm -o cson_F_hifi_phased.asm -t32 --h1 cson_f_hic_R1.fastq.gz --h2 cson_f_hic_R2.fastq.gz cson_f_hifi.filt.mitorm.fastq.gz
```