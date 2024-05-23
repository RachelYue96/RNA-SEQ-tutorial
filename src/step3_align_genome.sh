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
log_file="$wkdir/step3_align_genome.log"

commands=()
for sampleid in $samples
do 
	## 数据全路径 和 结果全路径 ##
	ddir="$wkdir/results/rm_rRNA/$sampleid"
	odir="$wkdir/results/hisat2_stringtie/$sampleid"

	## 创建结果目录 ##
	mkdir -p $odir

	## step 3. 比对genome ##
	cmd="echo \"$(date '+%Y-%m-%d %H:%M:%S') - Processing $sampleid\" >> $log_file; \
		hisat2 \
		--threads 16 \
		--dta \
		-x /share/data/reference/human/hg38/gene/UHRR/ref/genome_index/genome_GRCh38_110_index \
		-1 $ddir/${sampleid}.clean.fq.1.gz \
		-2 $ddir/${sampleid}.clean.fq.2.gz \
		--summary-file $odir/${sampleid}.tran_summary.txt \
		| samtools view \
			--threads 16 \
			--bam --no-PG /dev/stdin \
		| samtools sort \
			/dev/stdin \
			--threads 16 \
			-o $odir/${sampleid}.tran.sorted.bam \
		&& samtools index \
			$odir/${sampleid}.tran.sorted.bam \
		&& stringtie \
			$odir/${sampleid}.tran.sorted.bam \
			-p 16 \
			-G /share/data/reference/human/hg38/gene/UHRR/ref/Homo_sapiens.GRCh38.110.gtf \
			-e \
			-b $odir/Ballgown \
			-o $odir/${sampleid}.tran.gtf \
		>> $log_file 2>&1; \
        echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished processing $sampleid\" >> $log_file"

	commands+=("$cmd")
done

## 多线程分析全基因组比对
printf "%s\n" "${commands[@]}" | xargs -d '\n' -P $process_num -I {} bash -c '{}'

echo "$(date '+%Y-%m-%d %H:%M:%S') - All samples processed." >> $log_file
