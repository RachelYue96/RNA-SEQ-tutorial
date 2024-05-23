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
log_file="$wkdir/step4_count.log"

commands=()
for sampleid in $samples
do 
	## 数据全路径 和 结果全路径 ##
	ddir="$wkdir/results/hisat2_stringtie/$sampleid"

	## 创建结果目录 ##
	mkdir -p $ddir

    ## step 4. 基因和转录本水平的定量 ##
	cmd="echo \"$(date '+%Y-%m-%d %H:%M:%S') - Processing $sampleid\" >> $log_file; \
		featureCounts \
			-p \
			-T 8 \
			-f \
			-a /share/data/reference/human/hg38/gene/UHRR/ref/Homo_sapiens.GRCh38.110.mRNA.gtf \
			-o $ddir/${sampleid}_gene_level_counts \
			$ddir/${sampleid}.tran.sorted.bam && \
		featureCounts \
			-p \
			-T 8 \
			-g 'transcript_id' \
			-a /share/data/reference/human/hg38/gene/UHRR/ref/Homo_sapiens.GRCh38.110.mRNA.gtf \
			-o $ddir/${sampleid}_transcript_level_counts \
			$ddir/${sampleid}.tran.sorted.bam && \
		featureCounts \
			-p \
			-M \
			-O \
			-T 8 \
			-f \
			-a /share/data/reference/human/hg38/gene/UHRR/ref/Homo_sapiens.GRCh38.110.gtf \
			-o $ddir/${sampleid}_gene_ratio \
			$ddir/${sampleid}.tran.sorted.bam; \
		less $ddir/${sampleid}_transcript_level_counts | grep '^E' | awk '{print \$1\"\\t\"\$NF}' > $ddir/${sampleid}_transcript.counts.txt; \
		sed -i '1igene_id\tcounts' $ddir/${sampleid}_transcript.counts.txt; \
		less $ddir/${sampleid}_gene_level_counts | grep '^E' | awk '{print \$1\"\\t\"\$NF}' | awk '{sums[\$1] += \$2} END {for (i in sums) print i, sums[i]}' | sort -k1 > $ddir/${sampleid}_allgene.counts.txt; \
		sed -i '1igene_id\tcounts' $ddir/${sampleid}_allgene.counts.txt; \
		less $ddir/${sampleid}_gene_level_counts | grep '^E' | awk '{print \$1\"\\t\"\$(NF-1)}' > $ddir/${sampleid}_allgene.lengths.txt; \
		sed -i '1igene_id\tlengths' $ddir/${sampleid}_allgene.lengths.txt; \
		less $ddir/${sampleid}_transcript_level_counts | grep '^E' | cut -f1,6 > $ddir/${sampleid}_transcript.lengths.txt; \
		sed -i '1igene_id\tlengths' $ddir/${sampleid}_transcript.lengths.txt >> $log_file 2>&1; \
		echo \"$(date '+%Y-%m-%d %H:%M:%S') - Finished processing $sampleid\" >> $log_file"
    commands+=("$cmd")
done

## 多线程分析qc
printf "%s\n" "${commands[@]}" | xargs -d '\n' -P $process_num -I {} bash -c '{}'

mkdir -p $wkdir/results/RNA_counts
Rscript $wkdir/src/summary_counts.r

echo "$(date '+%Y-%m-%d %H:%M:%S') - All samples processed." >> $log_file


