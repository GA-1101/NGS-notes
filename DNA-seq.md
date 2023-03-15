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

### 数据预处理-碱基质量重校正

```
gatk BaseRecalibrator \
    -I YN20220865-M35_S1_bwa.sorted.markdup.bam \
    -R /home/zhanglei/ref/GCF_GRCh38.p14_genomic.fna \
    --known-sites /home/zhanglei/ref/VCF/1000G_phase1.snps.high_confidence.hg38.vcf \
    --known-sites /home/zhanglei/ref/VCF/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    -O YN20220865-M35_S1.recal_data.table
```

### Variants Calling

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