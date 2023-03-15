# DNA-seq

## Software prepare

### Anaconda
下载Anaconda包，并安装
conda用于提供生信分析的基础软件环境，包管理以及环境管理。

`https://www.anaconda.com/products/distribution#Downloads`

### Jupyter

`conda activate jupyter`

`jupyter notebook --ip=172.16.100.80 --port=8890 --no-browser`


## Pipeline

流程环境：
`conda activate DNA-seq`

### 质量控制
`fastqc /home/DATA/raw_data/*.fastq.gz -o /home/zhanglei/test/fastqc_out_dir/ -t 32`

### 建立参考基因组Index
`samtools faidx GCF_GRCh38.p14_genomic.fna`

`gatk CreateSequenceDictionary -R GCF_GRCh38.p14_genomic.fna`

### 序列比对
`bwa mem -t 32 /home/zhanglei/ref/Hg_38/GRCh38 ./raw_data/YN20220865-M35_S1_R1_001.fastq.gz ./raw_data/YN20220865-M35_S1_R2_001.fastq.gz > YN20220865-M35_S1_bwa.sam`

`samtools view -S -b -h YN20220865-M35_S1_bwa.sam -o YN20220865-M35_S1_bwa.bam`

### 数据预处理-排序

`samtools sort -@ 32 -O bam -o YN20220865-M35_S1_bwa.sorted.bam YN20220865-M35_S1_bwa.bam`

### 数据预处理-去重复

```
java -jar /home/zhanglei/software/picard.jar MarkDuplicates \
	I=YN20220865-M35_S1_bwa.sorted.bam \
	O=YN20220865-M35_S1_bwa.sorted.markdup.bam \
	M=YN20220865-M35_S1.markdup_metrics.txt
```

### 数据预处理-碱基质量重校正 - Base (Quality Score) Recalibration

```
gatk BaseRecalibrator \
    -I YN20220865-M35_S1_bwa.sorted.markdup.bam \
    -R /home/zhanglei/ref/GCF_GRCh38.p14_genomic.fna \
    --known-sites /home/zhanglei/ref/VCF/1000G_phase1.snps.high_confidence.hg38.vcf \
    --known-sites /home/zhanglei/ref/VCF/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    -O YN20220865-M35_S1.recal_data.table
```

### 变异检测 - Variants Calling

```
gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R /home/zhanglei/ref/GCF_GRCh38.p14_genomic.fna \
    -I fixed_YN20220865-M35_S1_bwa.sorted.markdup.bam \
    -A QualByDepth \
    -A RMSMappingQuality \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A Coverage \
    -O YN20220865-M35_S1.HC.vcf
```

### 变异检测 - Select Variants

SNP:
```
gatk SelectVariants \
    -R /home/zhanglei/ref/GCF_GRCh38.p14_genomic.fna \
    -V YN20220865-M35_S1.HC.vcf \
    --select-type-to-include SNP \
    -O YN20220865-M35_S1_snps.vcf
```

Indel:
```
gatk SelectVariants \
    -R /home/zhanglei/ref/GCF_GRCh38.p14_genomic.fna \
    -V YN20220865-M35_S1.HC.vcf \
    --select-type-to-include INDEL \
    -O YN20220865-M35_S1_indels.vcf
```

### 变异检测质控和过滤 - VQSR (Variant Quality Score Recalibration)

```
# SNPs VQSR
gatk VariantRecalibrator \
   -R <ucsc.hg3838.fasta> \
   -V <raw_snps.vcf> \
   --resource hapmap,known=false,training=true,truth=true,prior=15.0:hapmap_3.3.hg38.sites.vcf.gz \
   --resource omni,known=false,training=true,truth=false,prior=12.0:1000G_omni2.5.hg38.sites.vcf.gz \
   --resource 1000G,known=false,training=true,truth=false,prior=10.0:1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O <snps_recalibrate.recal> \
   --tranches-file <snps_recalibrate.tranches> \ 
   --rscript-file <snps_recalibrate.plots.R>
gatk ApplyVQSR \
   -R <ucsc.hg38.fasta> \
   -V <raw_snps.vcf> \
   -O <snps_recalibrate.vcf> \
   --truth-sensitivity-filter-level 99.5 \
   --tranches-file <snps_recalibrate.tranches> \
   --recal-file <snps_recalibrate.recal> \
   -mode SNP
##
# Indels VQSR
gatk VariantRecalibrator \
   -R <ucsc.hg38.fasta> \
   -V <raw_indels.vcf> \
   --maxGaussians 4 \
   -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf  \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_146.hg38.vcf\
   -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
   -mode INDEL \
   -O <indels_recalibrate.recal> \
   --tranches-file <indels_recalibrate.tranches> \
   --rscript-file <indels_recalibrate.plots.R>
gatk ApplyVQSR \
   -R <ucsc.hg38.fasta> \
   -V <raw_indels.vcf> \
   -O <indels_recalibrate.vcf> \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file <indels_recalibrate.tranches> \
   --recal-file <indels_recalibrate.recal> \
   -mode INDEL
```
