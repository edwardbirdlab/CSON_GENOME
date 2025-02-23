# CSON_GENOME
 This is a walk through of how I did the CSON Female Genome
 
 
 ## Adapter Removal
 
 apptainer version: 1.3.1 <br>
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

Minimap2 Version: 2.26 <br>
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

## Contamination Check with Blobtools

Minimap2 version: 2.26 <br>
Samtools Version: 1.17 <br>
Diamond Version: 2.0.15 <br>
Blobtools Version: 1.1 <br>
Uniprot & Blast DB Access Date: 2/21/2025

```
apptainer build blobtools_1.1.sif docker://ebird013/blobtools:1.1
apptainer build diamond_2.0.14.sif docker://bschiffthaler/diamond:2.0.14
apptainer exec ./minimap2_2.26.sif minimap2 -ax map-hifi cson_F_hifi_phased.asm.hic.hap1.p_ctg.fasta cson_f_hifi.filt.mitorm.fastq.gz > cson_f_coverage.sam
apptainer exec ./samtools_1.17.sif bash -c "samtools sort cson_f_coverage.sam -O bam -o cson_f_coverage.bam"
apptainer exec ./samtools_1.17.sif bash -c "samtools index -@ 16 cson_f_coverage.bam"
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz names.dmp nodes.dmp
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
apptainer exec ./diamond_2.0.14.sif diamond makedb --in uniprot_sprot.fasta -d uniprot_with_taxids --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp
apptainer exec ./diamond_2.0.14.sif diamond blastx --query cson_F_hifi_phased.asm.hic.hap1.p_ctg.fasta --db uniprot_with_taxids.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 16 > cson_uniprot_dblastx.out
mkdir blobout
apptainer exec ./blobtools_1.1.sif blobtools nodesdb --nodes nodes.dmp --names names.dmp
apptainer exec ./blobtools_1.1.sif blobtools create -i cson_F_hifi_phased.asm.hic.hap1.p_ctg.fasta -b cson_f_coverage.bam -t cson_uniprot_dblastx.out -o blobout/cson_f
apptainer exec ./blobtools_1.1.sif blobtools view -i blobout/cson_f.blobDB.json
apptainer exec ./blobtools_1.1.sif blobtools plot -i blobout/cson_f.blobDB.json
```

Insepct Blobtools output <br>

Blobtools Plot<br>

Make a list of the scaffolds you want to keep<br>

Blobtools Version: 1.1

```
apptainer exec ./blobtools_1.1.sif blobtools seqfilter -i cson_F_hifi_phased.asm.hic.hap1.p_ctg.fasta -l keep_scaffolds_blobtools.txt
```

Final Filtered Contigs at: cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fna <br>

## Juicer

Double Check enzyme site (DpnII)
```
head -n 500 Culicoides_sonorensis_F10_R1.fastq.gz | gunzip | grep -o 'GATCGATC' | wc -l
```

Juicer Version:
BWA Version: 
Samtools Version: 1.17

```
apptainer exec ./juicer_2.0.1.sif bwa index cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta
apptainer exec ./juicer_2.0.1.sif generate_site_positions.py DpnII modify_final cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta
samtools faidx cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta
aawk '/^>/ {if (seq) print name, length(seq); name = substr($1,2); seq=""; next} {seq = seq $1} END {if (seq) print name, length(seq)}' cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta > cson_F_chrom.size
mkdir fastq //Puts hi-c reads here
apptainer exec ./juicer_2.0.1.sif juicer.sh -g cson_F_hifi_phased -s DpnII -z cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta -y modify_final_DpnII.txt -p cson_F_chrom.sizes

wget https://raw.githubusercontent.com/aidenlab/3d-dna/master/utils/generate-assembly-file-from-fasta.awk
awk -f generate-assembly-file-from-fasta.awk cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta > cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.assembly
```

## HiCstuff

Double Check enzyme site (DpnII)
```
head -n 500 Culicoides_sonorensis_F10_R1.fastq.gz | gunzip | grep -o 'GATCGATC' | wc -l
```

HiCstuff version:

```
apptainer build hicstuff_3.2.4.sif docker://quay.io/biocontainers/hicstuff:3.2.4--pyhdfd78af_0
gunzip -f Culicoides_sonorensis_F10_R1.fastq.gz
gunzip -f Culicoides_sonorensis_F10_R2.fastq.gz
mkdir cson_f_hicstuff
apptainer exec ./hicstuff_3.2.4.sif hicstuff pipeline -t 16 -a minimap2 -e DpnII -o cson_f_hicstuff/ -g cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fna Culicoides_sonorensis_F10_R1.fastq Culicoides_sonorensis_F10_R2.fastq
awk '{print $1"\t"$2}' cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fna.fai > cson_f.chrom.sizes
```

## NF-Core hic

```
git clone https://github.com/nf-core/hic.git
cd hic/conf/
wget https://raw.githubusercontent.com/edwardbirdlab/rnaseq/master/conf/ceres.cfg
```
