# CSON GENOME Assembly & Annoations
 This is a walk through of how I did the CSON Female Genome <br>

 This is not intended as a turorial on how to do genome assembly, rather a log of what I did, follow at your own risk. <br>

 Between each major step I created a new clean directory, and symbolic linked relavent files to the new director. I reccomend this as some of these proceses generate a lot of files, and it will get unmanageable very quickly. <br>

 This work was completed on ceres on the scinet computing cluster, along side a workstation for certain steps
 
 
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
 
 MitoHiFi Version: 3.2.2 <br>
 Reference Mito Geome: BK065013.1
 
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

Depending on your sequencing it may be a good idea to run purge dupes at this point. In the the female cson genome our Hi-C data was not from the individaul that the HiFi data was from. While this lead to a nearly complete haplotype resolution a small amount of duplication was idntified during Hi-C scaffolding. Purge dupes may have helped with this.

## Contamination Check with Blobtools in Draft Assembly

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

## Create Draft Chromosome Assembly with 3d-dna

Double Check enzyme site (DpnII)
```
head -n 500 Culicoides_sonorensis_F10_R1.fastq.gz | gunzip | grep -o 'GATCGATC' | wc -l
```
<br>
Juicer Version: 2.0.1 <br>
Samtools Version: 1.17<br>
3d-dna Version: 190716

```
apptainer build juicer_2.0.1.sif
apptainer exec ./juicer_2.0.1.sif bwa index cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta
apptainer exec ./juicer_2.0.1.sif generate_site_positions.py DpnII modify_final cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta
apptainer exec ./samtools_1.17.sif samtools faidx cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta
awk '/^>/ {if (seq) print name, length(seq); name = substr($1,2); seq=""; next} {seq = seq $1} END {if (seq) print name, length(seq)}' cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta > cson_F_chrom.size
wget https://raw.githubusercontent.com/aidenlab/3d-dna/master/utils/generate-assembly-file-from-fasta.awk
awk -f generate-assembly-file-from-fasta.awk cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta > cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.assembly
mkdir fastq #Put hi-c reads in this folder. Make sure they have the correct naming scheme *_R1.fastq.gz
apptainer exec ./juicer_2.0.1.sif juicer.sh -g cson_F_hifi_phased -s DpnII -z cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta -y modify_final_DpnII.txt -p cson_F_chrom.sizes --assembly cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.assembly
mkdir 3dddna
cd 3ddna
apptainer build 3d_dna.sif quay.io/biocontainers/3d-dna:201008--hdfd78af_0
apptainer exec ./3d_dna.sif 3d-dna --assembly ../cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.assembly ../aligned/merged_nodups.txt
```
## Manual Curation of Draft Chromosome Assemlby in Juicebox

Juicebox Version: 2.15
3d-dna Version: 190716

For this step we will download 2 of the output files for 3d-dna and load them into juicebox. Files cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.rawchromo...assembly and cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.rawchromo...hic.<br>

These files were loaded into juice box and manually scaffolded according to the instrdtions in the Genome Assembly Cookbook from (aidenlab) - (https://aidenlab.org/assembly/manual_180322.pdf) <br>

Final HiC Map: <br>

Relatively minor changes were made to the map that was ouput from 3d-dna. Our map shows clear chromosome level assembly, and has relatively few unplaced scaffolds. A couple of these unplaced scaffolds show contact with some chromosome, but is ambigious, and were noth placed. Many of the smaller debris peices (pieces cut out from missassembly) show no contact points, even to themselves.  <br>

## Producing Final Gapped Chromosome Assembly

The edited assembly map was exported from jucebox and upload. The map was edited on a windows compter, so I convert line endings just to be safe. <br>

Finally we will output our final gapped chromosome by running 3d-dna one last time
```
cat cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.rawchrom.review.assembly | sed 's/\r$//' > cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.rawchrom.fix.review.assembly
apptainer exec /usr/local/share/3d-dna/run-asm-pipeline-post-review.sh -r cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.rawchrom.fix.review.assembly cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered.fasta merged_nodups.txt
```
We can also inspect the final Hi-C graph to make sure it looks as we expect <br>

## NCBI - Foreign Containiment Screen on Final Chromosome Assembly

The containiment check eariler was a quick and dirty protein blast using the uniprot database. This can lead to some odd classifications, and this time we want to be absolutely sure if we should remove sesquences or not.<br>

FCS version: 3.0.4 <br>
FCS-gx Datase Version: build:Jan 19 2023 15:50:20; git:v0.3.0-151-g9aad15db - Accessed on 2/26/2025 <br>

```
curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-gx.sif -Lo fcs-gx.sif
LOCAL_DB="./fcs_gx_db"
curl -LO https://github.com/peak/s5cmd/releases/download/v2.0.0/s5cmd_2.0.0_Linux-64bit.tar.gz
tar -xvf s5cmd_2.0.0_Linux-64bit.tar.gz
./s5cmd  --no-sign-request cp  --part-size 50  --concurrency 50 s3://ncbi-fcs-gx/gxdb/latest/all.* $LOCAL_DB
apptainer exec ./fcs-gx.sif python3 /app/bin/run_gx --fasta cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered_HiC.fasta --out-dir ./gx_out --gx-db ./fcs_gx_db --tax-id 179676
```

There was no containimation detected by FCS in out final genome. However, there were sevaral debris pieces that had no identifiable hit in gx_dx. This included HiC_scaffold_6, HiC_scaffold_13, HiC_scaffold_14, HiC_scaffold_17, HiC_scaffold_18, HiC_scaffold_20, and HiC_scaffold_28.

## Checking for duplication in debris

```
awk '/^>/ {count++} count<=3' cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered_HiC.fasta > cson_F_chromosomes.fasta
awk '/^>/ {count++} count>3' cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered_HiC.fasta > cson_F_debris.fasta
apptainer exec ./mummer4.sif dnadiff cson_F_chromosomes.fasta cson_F_debris.fasta
```
## Rules for Removing Unplaced Scaffolds

Debris sequences were removed only if they met one (or more) of the following conditions:

### Type I Removal: Haplotype contamination - observed by Hi-C map
 - This type of debis occured via haplotype contamination. In this type the deduplicatied HiC data mapped better to the unplaced debris that the scaffolded chromosome. This looks like a lack of hic evidend in the chromosome, but a stong, linear interaction between part/all of the debris sample. Below is a example of what this looks like on the HiC map. This was considered enough evidence to remove these scaffolds.
   _Insert example here_
### Type II Removal: Haplotype contamination - observed by Hi-C & Mummer
 - This type of haplotype contaomination is harder to determine than the above, so more evidence has been incorperated. This is the reverse, where the unplaced scaffold has very little to no HiC contacts. On possibile reson for this is that it is a haplotig, and the HiC data mapped better to the chromosome integrated sequence that this unplaced scaffold. However, since the HiC data was from the mother, and not the HiFi individual we wanted additional evidence to back this up. We ran Mummer dna diff, and looked for sequences that had a very high coverage query coverage compared to the chromosomes (Greater than 90%).
### Type III Removal: Small with No Annotations
 - HiC scaffolding left several small, unplaced scaffolds due either manual, or 3ddna miscaffolding snipping. These are very small, usually 1kb, or 5kb. For this reason we will remove any sequenced less than 10kb that had no predicted annotation from eGAPx.
### Type IV Removal: Unidentified sequence with No Annotations
 - While Foreign Containiment Screen did not find any evidence of contamination in our sequences, it was unable to classify some of our sequences. Therefore any sequences that FCS was not able to identify _and_ get no annoations from eGAPx will be removed.
### Type V Removal:

## Removal of Unplace Scaffolds Phase 1 (Before prelim eGAPx annotation)

Several Scaffolds met one of the conditions before running eGAPx (Phase 1 recommendation in spread sheet) were removed. <br>
 This included: <br>
 HiC_scaffold_10 - Type I<br>
 HiC_scaffold_31 - Type II<br>
 HiC_scaffold_32 - Type I & II



## Renamming Scaffolds

Chromosome were nammed in order of size:<br>
chr1, chr2, chr3<br>

Unplaced scaffolds that had Hi-C contact evidence with a certain chromosome (ambiguous on exact position) were labeled (number by size):<br>
chr#_unplaced\_#<br>

Unplaced scaffolds that had no clear Hi-C contatct evidnece with a _single_ chromosome were labeled (number by size):<br>
unplaced_#<br>

The fasta was nammed cson_F_hifi_phased.asm.hic.hap1.p_ctg.filtered_HiC.phase1.fasta after removal of the 3 scaffolds and name change

## Gathering SRA RNAseq Data

Now we will gather and trim rnaseq data from SRA, and Trim data, and inspect for QC. All accessions for tasid 179676<br>

After getting a the SRA accession sheet from NCBi we will use edwardbirdlab/SRA-FETCH

```
nextflow run -resume /project/culicoides/bird_projects/pipelines/SRA-FETCH -c /project/culicoides/bird_projects/pipelines/SRA-FETCH/configs/ceres/ceres.cfg \
    --workflow_opt paired \
    --project_name cson_rnaseq \
    --sample_sheet ./cson_rnaseq_sra_acc.csv \
    --apikey e10498aa1ed5c5feb58f5ddcec0cfe381a09 \
    --fastp_minlen 50 \
    -N edwardbird@ksu.edu

nextflow run -resume /project/culicoides/bird_projects/pipelines/SRA-FETCH -c /project/culicoides/bird_projects/pipelines/SRA-FETCH/configs/ceres/ceres.cfg \
    --workflow_opt multiqc \
    --project_name cson_rnaseq \
    --sample_sheet ./cson_rnaseq_sra_acc.csv \
    --apikey e10498aa1ed5c5feb58f5ddcec0cfe381a09 \
    --fastp_minlen 50 \
    -N edwardbird@ksu.edu
```
Investigate MultiQC to determine if RNASeq data should be used for genome annotaion <br>

All samples except (insert sample here) were determined to be of sufficent quality.

## Preping Life History RNAseq Data

Now we will also trim the MLH data and ispect its QC

I used the Shortread QC module from edwardbirdlab/BALROG-MSR

```
nextflow run -resume /90daydata/shared/edwardbird/Pipelines/BALROG-MSR -c /90daydata/shared/edwardbird/Pipelines/BALROG-MSR/configs/ceres/ceres.cfg \
    --workflow_opt sr_qc_only \
    --project_name mlh_trim \
    --sample_sheet ./mlh_samplesheet.csv \
    --fastp_minlen 50 \
    -N edwardbird@ksu.edu

nextflow run -resume /90daydata/shared/edwardbird/Pipelines/BALROG-MSR -c /90daydata/shared/edwardbird/Pipelines/BALROG-MSR/configs/ceres/ceres.cfg \
    --workflow_opt multiqc \
    --project_name mlh_trim \
    --sample_sheet ./mlh_samplesheet.csv \
    --fastp_minlen 50 \
    -N edwardbird@ksu.edu
```
Inspected QC and determine all samples were of sufficent quality.

## Running eGAPx Genome Annotation (Prelim Annotaiton)

Scinet was blocked from the NCBI FTP API at this time, so I had to donload the local cache and upload it to run eGAPx in offline mode <br>

```
source /project/culicoides/bird_projects/conda/epagx/bin/activate
python3 ./ui/egapx.py /90daydata/shared/edwardbird/genome_reassembly/egapx/egapx/input_cson_all.yaml -e slurm -o cson_all -lc ../local_cache
```

## Removal of Unplace Scaffolds Phase 12 (After prelim eGAPx annotation)

## Running eGAPx Genom Annotaion (Final Assembely)

## Running Earl Grey
