# Alignment with alternative reference genomes

In this tutorial, we will demonstrate an alignment pipeline that utilized an alternative reference genome and `levioSAM` to enchance accuracy.

## Software and datasets involved in this pipeline
Software:
- levioSAM
- mason2 simulator
- bcftools
- samtools
- Bowtie 2 or Bwa-mem

Datasets:
- GRCh38
- Small variant calls from the 1000 Genomes Project
- A chromosome name map file (map chromosome name "Z" to "chrZ")

```shell
DIR_DATA="resources"
# GRCh38
wget -P ${DIR_DATA}/ ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
bgzip -d ${DIR_DATA}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx ${DIR_DATA}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna chr21 > ${DIR_DATA}/chr21.fa
# chromosome name map file
wget -P ${DIR_DATA}/ https://raw.githubusercontent.com/langmead-lab/reference_flow/master/resources/GRCh38.chrom_map
# 1KGP genotypes
wget -P ${DIR_DATA}/1kg_vcf/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
```



## Simulate data

We will use a simulated dataset in this tutorial. The advantage of simulated sequencing data is the advent of a ground truth.

We use the genotypes from the first haplotype of individual NA12878 from the 1000 Genomes Project (1KGP). We only use chr21 data for simplicity.

```sh
SAMPLE="NA12878"
bcftools view -s ${SAMPLE} ${DIR_DATA}/1kg_vcf/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | bcftools annotate --rename-chrs ${DIR_DATA}/GRCh38.chrom_map -O z -o chr21-${SAMPLE}.vcf.gz; bcftools index chr21-${SAMPLE}.vcf.gz
# Generate NA12878 chr21 hapA. 
# This will be used as the personalized reference genome later.
bcftools consensus -f ${DIR_DATA}/chr21.fa -H 1 -o ${SAMPLE}-chr21-hapA.fa chr21-${SAMPLE}.vcf.gz
# Simulate reads using mason2
mason_simulator -ir ${SAMPLE}-chr21-hapA.fa -o ${SAMPLE}-chr21-hapA.1.fq -or ${SAMPLE}-chr21-hapA.2.fq -oa ${SAMPLE}-chr21-hapA-mason.sam -n 100000 --num-threads 8
```



## Build major-allele reference

```sh
bcftools view -O z -q 0.5000001 -G ${DIR_DATA}/1kg_vcf/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | bcftools annotate --rename-chrs ${DIR_DATA}/GRCh38.chrom_map -O z -o chr21-major.vcf.gz 
bcftools index chr21-major.vcf.gz
bcftools consensus -f ${DIR_DATA}/chr21.fa -o chr21-major.fa chr21-major.vcf.gz
```



## Build indexes for Bowtie 2/bwa-mem

```sh
# Bowtie 2
mkdir bt2_indexes
bowtie2-build --threads 8 chr21-major.fa bt2_indexes/chr21-major.fa
bowtie2-build --threads 8 ${DIR_DATA}/chr21.fa bt2_indexes/chr21.fa
bowtie2-build --threads 8 ${SAMPLE}-chr21-hapA.fa bt2_indexes/${SAMPLE}-chr21-hapA.fa
# bwa mem
mkdir bwa_indexes
bwa index -p bwa_indexes/chr21-major.fa chr21-major.fa
bwa index -p bwa_indexes/chr21.fa ${DIR_DATA}/chr21.fa
bwa index -p bwa_indexes/${SAMPLE}-chr21-hapA.fa ${SAMPLE}-chr21-hapA.fa
```


