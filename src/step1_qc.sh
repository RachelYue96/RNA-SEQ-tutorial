#!/bin/bash

## 参数 ##
## wkdir - working directory
## process_num - 同时运行的进程数
wkdir="/share/result/sequencer/salus/video_example/RNA-seq"
process_num=3
cd $wkdir

## 样本名 ##
samples=$(find ./data -name '*_R1.fastq.gz' | sed 's/_R1.fastq.gz//' | xargs -n 1 basename)

## 数据全路径 和 结果全路径 ##
ddir="$wkdir/data"

## log文件
log_file="$wkdir/step1_qc.log"

commands=()
for sampleid in $samples
do 
    odir="$wkdir/results/fastp/$sampleid"
    mkdir -p $odir  # 创建结果目录

    cmd="echo \"$(date '+%Y-%m-%d %H:%M:%S') - Processing $sampleid\" >> $log_file; \
         fastp \
         --in1 $ddir/${sampleid}_R1.fastq.gz \
         --out1 $odir/${sampleid}_1.clean.fastq.gz \
         --in2 $ddir/${sampleid}_R2.fastq.gz \
         --out2 $odir/${sampleid}_2.clean.fastq.gz \
         --adapter_fasta /share/data/reference/human/ALK/alk_adapter.fa \
         --detect_adapter_for_pe \
         --qualified_quality_phred 20 \
         --unqualified_percent_limit 20 \
         --length_required 150 \
         --n_base_limit 5 \
         --low_complexity_filter \
         --complexity_threshold 30 \
         --json $odir/${sampleid}.fastp.json \
         --html $odir/${sampleid}.fastp.html \
         >> $log_file 2>&1; \
         echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished processing $sampleid\" >> $log_file"
    
    commands+=("$cmd")
done

## 多线程分析qc
printf "%s\n" "${commands[@]}" | xargs -d '\n' -P $process_num -I {} bash -c '{}'

echo "$(date '+%Y-%m-%d %H:%M:%S') - All samples processed." >> $log_file