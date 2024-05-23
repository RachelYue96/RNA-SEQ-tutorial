#!/bin/bash

## 参数 ##
## wkdir - working directory
## process_num - 同时运行的进程数
wkdir="/share/result/sequencer/salus/video_example/RNA-seq"
process_num=3
cd $wkdir

## 样本名 ##
samples=$(find ./data -name '*_R1.fastq.gz' | sed 's/_R1.fastq.gz//' | xargs -n 1 basename)

# log文件
log_file="$wkdir/step2_rm_rRNA.log"

commands=()
for sampleid in $samples
do 
	## 数据全路径 和 结果全路径 ##
	ddir="$wkdir/results/fastp/$sampleid"
	odir="$wkdir/results/rm_rRNA/$sampleid"

	## 创建结果目录 ##
	mkdir -p $odir

	## step 2. 去rRNA ##
	cmd="echo \"$(date '+%Y-%m-%d %H:%M:%S') - Processing $sampleid\" >> $log_file; \
		bowtie2 \
		--very-sensitive-local \
		--no-unal \
		--threads 16 \
		-x /share/data/reference/human/grch38/rRNA/GRCh38.p14_rRNA_rm16dup.fasta \
		-1 $ddir/${sampleid}_1.clean.fastq.gz \
		-2 $ddir/${sampleid}_2.clean.fastq.gz \
		--un-conc-gz $odir/${sampleid}.clean.fq.gz 2> $odir/${sampleid}.map2rRNAstats.txt \
		| samtools view \
			--no-PG \
			--bam \
			--output $odir/${sampleid}.map2rRNA.bam - \
		>> $log_file 2>&1; \
        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished processing $sampleid\" >> $log_file"
	
	commands+=("$cmd")
done

## 多线程分析去rRNA
printf "%s\n" "${commands[@]}" | xargs -d '\n' -P $process_num -I {} bash -c '{}'

echo "$(date '+%Y-%m-%d %H:%M:%S') - All samples processed." >> $log_file
