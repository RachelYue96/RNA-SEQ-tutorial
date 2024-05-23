# RNA测序数据分析流程 - 教程详解
在这个RNA测序数据分析的教程中，我们将以六组PE100的测序FASTQ数据（样本类型：6例人宫颈癌培养细胞样本）为例介绍RNA测序数据的差异基因表达常规分析流程。

其中，FASTQ数据来源：采用特有的芯片结合可逆末端终止碱基的测序技术，通过光学系统对 DNA 分子群的荧光信号进行图像采集，得到的图像信息通过 Basecall转换为待测碱基的序列信息。

## 分析全流程 ##

```bash
## 第一步：数据质控 ##
bash src/step1_qc.sh

## 第二步：去rRNA ##
bash src/step2_rm_rRNA.sh

## 第三步：比对到参考基因组
bash src/step3_align_genome.sh

## 第四步：基因和转路本定量
bash src/step4_count.sh

## 第五步：画图，差异基因分析及通路分析
bash src/step5_plots.sh
```

## 使用到的开源软件 ##
- [FASTP](https://github.com/OpenGene/fastp)

- [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- [HISAT2](https://daehwankimlab.github.io/hisat2/)

- [samtools](https://github.com/samtools/samtools)

- [StringTie](https://ccb.jhu.edu/software/stringtie/)

- [featureCounts](https://subread.sourceforge.net/featureCounts.html)

- [R](https://www.r-project.org/)

## 第一步：数据质控 ##
数据质量控制是RNA-Seq分析的关键步骤之一，这里使用软件fastp对原始测序FASTQ数据的测序reads进行质量评估和过滤，以确保我们分析的FASTQ数据质量良好，减少后续分析中的误差。

*完整代码参考 [github](https://github.com/RachelYue96/RNA-SEQ-tutorial) `src/step1_qc.sh` 脚本*

### 核心代码 ###
```bash
fastp \
    --in1 $ddir/${sampleid}_R1.fastq.gz \ ## 双端R1输入文件
    --out1 $odir/${sampleid}_1.clean.fastq.gz \ ## 双端R1输出文件
    --in2 $ddir/${sampleid}_R2.fastq.gz \ ## 双端R2输入文件
    --out2 $odir/${sampleid}_2.clean.fastq.gz \ ## 双端R2输出文件
    --adapter_fasta /share/data/reference/human/ALK/alk_adapter.fa \
    --detect_adapter_for_pe \ ## 自动检测双端数据接头
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 20 \
    --length_required 150 \
    --n_base_limit 5 \ ## 设置一个序列中允许的最大 N 碱基数
    --low_complexity_filter \
    --complexity_threshold 30 \
    --json $odir/${sampleid}.fastp.json \ ## json结果文件
    --html $odir/${sampleid}.fastp.html \ ## html结果文件
    >> $log_file 2>&1 ## 导出日志
```

### 结果文件 ###
- 过滤后的fastq双端数据
- html结果文件
- json结果文件 </br>

![Page Preview](images/json_result.png)

## 第二步：去rRNA ##
rRNA在RNA-Seq分析中通常不是我们关心的目标，这里使用软件bowtie2将测序FASTQ数据与已知的rRNA序列进行比对，再过滤掉比对上的reads

*完整代码参考 [github](https://github.com/RachelYue96/RNA-SEQ-tutorial) `src/step2_rm_rRNA.sh` 脚本*

### 下载rRNA参考基因组 ###
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz

# extract rRNA seqs
zcat GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz |grep "^>" |grep "gbkey=rRNA" | awk '{print $1}'|sed 's/>//g' > rRNA_id.list

seqkit grep -f rRNA_id.list GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz > GRCh38.p14_rRNA.fasta #~3s

# 前17条 5S rRNA序列除了第9条有一个碱基之外，其他序列是一模一样的，去冗余只保留一条
le rRNA_id.list |sed '2,17d' > rRNA_rm16dup.id.list
seqkit grep -f rRNA_rm16dup.id.list GRCh38.p14_rRNA.fasta > GRCh38.p14_rRNA_rm16dup.fasta

# fai and dict
samtools faidx GRCh38.p14_rRNA_rm16dup.fasta
samtools dict GRCh38.p14_rRNA_rm16dup.fasta > GRCh38.p14_rRNA_rm16dup.dict

# bowtie2 index 生成索引文件
bowtie2-build GRCh38.p14_rRNA_rm16dup.fasta GRCh38.p14_rRNA_rm16dup.fasta #~2s
```
### 核心代码 ###
```bash
bowtie2 \
    --very-sensitive-local \
    --no-unal \
    --threads 16 \ ## 线程数
    -x /share/data/reference/human/grch38/rRNA/GRCh38.p14_rRNA_rm16dup.fasta \ ## 参考基因组
    -1 $ddir/${sampleid}_1.clean.fastq.gz \ ## fastp输出的R1结果文件
    -2 $ddir/${sampleid}_2.clean.fastq.gz \ ## fastp输出的R2结果文件
    --un-conc-gz $odir/${sampleid}.clean.fq.gz 2> $odir/${sampleid}.map2rRNAstats.txt \
    | samtools view \
            --no-PG \
            --bam \ ## 输出bam格式的结果文件
            --output $odir/${sampleid}.map2rRNA.bam - \ ## 设定结果文件名
    >> $log_file 2>&1 ## 导出日志文件
```





