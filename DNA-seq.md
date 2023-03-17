# DNA-seq

## Software prepare

### Anaconda
下载Anaconda包，并安装
conda用于提供生信分析的基础软件环境，包管理以及环境管理。

`https://www.anaconda.com/products/distribution#Downloads`

### Jupyter

安装Jupyter Notebook，Jupyter Notebook是用于数据分析处理工作的一种集成开发环境，可以根据不同的conda环境，启动相应的内核，用于编写和调试脚本。

`conda activate jupyter`

`jupyter notebook --ip=<ip> --port=<port> --no-browser`


## Pipeline

流程环境：
`conda activate <environment>`

### 质量控制

在做read质量值分析的时候，FastQC并不单独查看具体某一条read中碱基的质量值，而是将Fastq文件中所有的read数据都综合起来一起分析。
`fastqc <fastq path> -o <output directory> -t <tread number>`

### 建立参考基因组Index


`samtools faidx <ref.fasta>`

`gatk CreateSequenceDictionary -R <ref.fasta>`

### 序列比对


`bwa mem -t 32 <ref files path> <Sample_Read1.fastq.gz> <Sample_Read2.fastq.gz> > Sample_bwa.sam`

`samtools view -S -b -h Sample_bwa.sam -o Sample_bwa.bam`

### 数据预处理-排序

`samtools sort -@ 32 -O bam -o Sample_bwa.sorted.bam Sample_bwa.bam`

### 数据预处理-去重复

```
java -jar picard.jar MarkDuplicates \
	I=Sample_bwa.sorted.bam \
	O=Sample_bwa.sorted.markdup.bam \
	M=Sample.markdup_metrics.txt
```

### 数据预处理-添加ReadGroups（可选）

某些数据可能因为bam文件中缺少ReadGroups，后续变异检测步骤报错： `the sample list cannot be null or empty`
此时可以通过picard的 `AddOrReplaceReadGroups` 功能手动添加这一信息。

相关解决方案参考自[GATK论坛](https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information)。

```
java -jar picard.jar AddOrReplaceReadGroups \
    I=Sample_bwa.sorted.markdup.bam \
    O=fixed_Sample_bwa.sorted.markdup.bam \
    SORT_ORDER=coordinate \
    RGID=foo \
    RGLB=bar \
    RGPL=illumina \
	RGPU=unit1 \
    RGSM=Sample1 \
    CREATE_INDEX=True
```


### 数据预处理-碱基质量重校正 - Base (Quality Score) Recalibration

此步骤可能非必要，请参照后续更新内容。

```
gatk BaseRecalibrator \
    -I Sample_bwa.sorted.markdup.bam \
    -R /home/zhanglei/ref/GCF_GRCh38.p14_genomic.fna \
    --known-sites /home/zhanglei/ref/VCF/1000G_phase1.snps.high_confidence.hg38.vcf \
    --known-sites /home/zhanglei/ref/VCF/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    -O Sample.recal_data.table
```

### 变异检测 - Variants Calling

```
gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R /home/zhanglei/ref/GCF_GRCh38.p14_genomic.fna \
    -I fixed_Sample_bwa.sorted.markdup.bam \
    -A QualByDepth \
    -A RMSMappingQuality \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A Coverage \
    -O Sample.HC.vcf
```

### 变异检测 - Select Variants

SNP:
```
gatk SelectVariants \
    -R <ref.fasta> \
    -V Sample.HC.vcf \
    --select-type-to-include SNP \
    -O Sample_snps.vcf
```

Indel:
```
gatk SelectVariants \
    -R <ref.fasta> \
    -V Sample.HC.vcf \
    --select-type-to-include INDEL \
    -O Sample_indels.vcf
```

### 变异检测质控和过滤 - VQSR (Variant Quality Score Recalibration)

```
# SNPs VQSR
gatk VariantRecalibrator \
   -R <ucsc.hg38.fasta> \
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
