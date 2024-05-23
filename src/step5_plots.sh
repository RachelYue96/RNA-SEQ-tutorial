#!/bin/bash

## 参数 ##
## wkdir - working directory
## process_num - 同时运行的进程数
wkdir="/share/result/sequencer/salus/video_example/RNA-seq"
bed_file="/share/data/reference/human/hg38/gene/UHRR/ref/Homo_sapiens.GRCh38.110.bed"
process_num=3
cd $wkdir

## 样本名 ##
samples=$(find ./data -name '*_R1.fastq.gz' | sed 's/_R1.fastq.gz//' | xargs -n 1 basename)

# log文件
log_file="$wkdir/step5_plots.log"

commands=()
for sampleid in $samples
do
    ## 数据全路径 和 结果全路径 ##
    odir="$wkdir/results/hisat2_stringtie/$sampleid"

    cmd="echo \"$(date '+%Y-%m-%d %H:%M:%S') - Processing $sampleid\" >> $log_file; \
        samtools view -b -L $bed_file \
        $odir/$sampleid.tran.sorted.bam > $odir/$sampleid.tran.extracted.bam && \
        samtools sort $odir/$sampleid.tran.extracted.bam -o $odir/$sampleid.tran.extracted.sorted.bam && \
        samtools index $odir/$sampleid.tran.extracted.sorted.bam && \
        geneBody_coverage.py -r $bed_file \
        -i $odir/$sampleid.tran.extracted.sorted.bam -o $odir/output -l 150 -f pdf && \
        rm $odir/$sampleid.tran.extracted.bam >> $log_file 2>&1; \
		echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished processing $sampleid\" >> $log_file"

    commands+=("$cmd")
done

## 多线程分析qc
printf "%s\n" "${commands[@]}" | xargs -d '\n' -P $process_num -I {} bash -c '{}'

## 找差异基因，画PCA图
Rscript $wkdir/src/general_plots.r

## 通路分析
Rscript $wkdir/src/pathway_analysis.r

echo "$(date '+%Y-%m-%d %H:%M:%S') - All samples processed." >> $log_file