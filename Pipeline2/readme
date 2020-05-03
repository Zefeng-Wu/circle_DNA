#!/bin/bash

#bash readme.sh
date

input_dir=/home/wuzefeng/Myresearch/ZZ/2020.3.28
work_dir=/home/wuzefeng/Myresearch/ZZ/2020.3.28
ref_genome=~/Myresearch/Genomes/DNA/GRCh38.primary_assembly.genome.fa
#faidx -x $ref_genome/GRCh38.primary_assembly.genome.fa # split chromsomes
ref_chromosomes=/home/wuzefeng/Myresearch/Genomes/DNA/human_chromsomes

mkdir $work_dir/1repeated_reads
mkdir $work_dir/1non-repeated_reads

# Tidehunter
for m in $(ls $input_dir/*.fq);do echo $m; filename=$(basename $m); TideHunter -o $work_dir/1repeated_reads/1_${filename%.fq}.cons.fa -t 3 $input_dir/$filename;done

## output fastq files with repeated reads and non-repeated reads
##################################################################

#grep ">" 1_*.cons.fa | awk '{u=split($1,m,"_");print m[1]}' | sort| uniq >repeated_reads/repeated_reads_id.txt
#cat repeated_reads/repeated_reads_id.txt | while read LINE; do echo $LINE; cat *.fq | awk 'BEGIN{RS="@";FS="\n"}; $1~/'$LINE'/{print ">"$1"\n"$2; exit}' >> repeated_reads.fq ;done
for m in $(ls $input_dir/*.fq);do echo $m; filename=$(basename $m); Rscript Fq.extract.R $filename $work_dir/1repeated_reads/1_${filename%.fq}.cons.fa; done
###################################################################

cd $work_dir/1repeated_reads
# minimap2
#for m in $(ls 1_*.fa);do echo $m; minimap2 -ax map-ont  $ref_genome $m -t 1 --sam-hit-only --secondary=no -o 2_${m:2:-3}.minimap2.sam;done

for m in $(ls 1_*.cons.fa) 
do 
	for n in $(ls $ref_chromosomes/*.fa)
	do (echo $m,$n
		prefix=$(basename $n)
    	minimap2 -ax map-ont  $n $m -t 1 --sam-hit-only --secondary=no -o 2_${m:2:-3}.$prefix.minimap2.sam
       )
    done
done

## get sample names
samples_names=$(ls *.sam | awk 'BEGIN{FS="."}{print$1}'| sed 's/2_//g'| sort| uniq)

## merge result from same sample
for sample in $samples_names; do cat *$sample*.sam > 3_$sample.sam;done

## remove temp results
rm 2_*.sam

## shorten reads name to meet samtools sort
for m in $(ls 3_*.sam);do awk 'OFS="\t"{if($1~/^@/)print $0;  else {u=split($1,m,"_"); $1=m[1]"_"m[2];  print $0}}' $m >4_${m:2:-4}.modified.sam;done

## reads count for each chromsomes
for m in $(ls *.modified.sam);do grep -v @ $m | cut -f1,3 | awk '{u=split($1,m,"_");print m[1], $2}' | sort | uniq| cut -d " " -f2 | sort | uniq -c | sort -nr > reads2chr.txt;done

## header modify and sort sam
for sam in $(ls 4_*.sam); do sed '/^@PG/d' $sam | sort -k1r| samtools view -S -b |samtools sort -@ 3 > 5_${sam:2:-4}.sorted.bam;done

## bam index
for m in $(ls 5_*.bam);do samtools index $m;done


## output fastq files with mapped consensus reads
##################################################################
mkdir mapping
for m in $(ls $input_dir/*.fq);do echo $m; filename=$(basename $m); Rscript ../Fq.extract2.R ../$filename 1_${filename%.fq}.cons.fa 5_${filename%.fq}.modified.sorted.bam; done
###################################################################



## bam to bedgraph
for m in $(ls 5_*.bam);do bamCoverage -b $m --outFileFormat bedgraph -o 6_${m:2:-4}.bedgraph -p 3 --binSize 50; done

## merge bam to get uniq region

for m in $(ls 5_*.sorted.bam); do bedtools merge -i $m > 7_${m:2:-4}.merge_bam.bed;done

for m in $(ls 7_*.bed);do bedtools coverage -a $m -b 5_${m:2:-14}.bam | awk '$4>=1 {print $0}'> 8_${m:2:-4}.reads_count.txt; done #每个bed区间的reads数和总碱基数

date

###################################################################
###################################################################
##  simialr analysis for non-repeated reads 
###################################################################
###################################################################

cd $work_dir/1non-repeated_reads
# minimap2
#for m in $(ls 1_*.fa);do echo $m; minimap2 -ax map-ont  $ref_genome $m -t 1 --sam-hit-only --secondary=no -o 2_${m:2:-3}.minimap2.sam;done

for m in $(ls 1_*.fa) 
do 
	for n in $(ls $ref_chromosomes/*.fa)
	do (echo $m,$n
		prefix=$(basename $n)
    	minimap2 -ax map-ont  $n $m -t 1 --sam-hit-only --secondary=no -o 2_${m:2:-3}.$prefix.minimap2.sam
       )
    done
done

## get sample names
samples_names=$(ls *.sam | awk 'BEGIN{FS="."}{print$1}'| sed 's/2_//g'| sort| uniq)

## merge result from same sample
for sample in $samples_names; do cat *$sample*.sam > 3_$sample.sam;done

## remove temp results
rm 2_*.sam

## shorten reads name to meet samtools sort
for m in $(ls 3_*.sam);do awk 'OFS="\t"{if($1~/^@/)print $0;  else {u=split($1,m,"_"); $1=m[1]"_"m[2];  print $0}}' $m >4_${m:2:-4}.modified.sam;done

## reads count for each chromsomes
for m in $(ls *.modified.sam);do grep -v @ $m | cut -f1,3 | awk '{u=split($1,m,"_");print m[1], $2}' | sort | uniq| cut -d " " -f2 | sort | uniq -c | sort -nr > reads2chr.txt;done

## header modify and sort sam
for sam in $(ls 4_*.sam); do sed '/^@PG/d' $sam | sort -k1r| samtools view -S -b |samtools sort -@ 3 > 5_${sam:2:-4}.sorted.bam;done

## bam index
for m in $(ls 5_*.bam);do samtools index $m;done


## output fastq files with mapped unrepeated reads
##################################################################
mkdir mapping
for m in $(ls $input_dir/*.fq);do echo $m; filename=$(basename $m); Rscript ../Fq.extract2.R ../$filename 1_${filename%.fq}.original_reads.fa 5_${filename%.fq}.modified.sorted.bam; done
###################################################################

## bam to bedgraph
for m in $(ls 5_*.bam);do bamCoverage -b $m --outFileFormat bedgraph -o 6_${m:2:-4}.bedgraph -p 3 --binSize 50; done

## merge bam to get uniq region

for m in $(ls 5_*.sorted.bam); do bedtools merge -i $m > 7_${m:2:-4}.merge_bam.bed;done

for m in $(ls 7_*.bed);do bedtools coverage -a $m -b 5_${m:2:-14}.bam | awk '$4>=1 {print $0}'> 8_${m:2:-4}.reads_count.txt; done #每个bed区间的reads数和总碱基数



#################################################################
Reads to transposons analysis
#################################################################
mkdir $work_dir/2TE_reads
mkdir $work_dir/2non_TE_reads

TE_seqs=/home/wuzefeng/Myresearch/ZZ/hg19.repBase.fa

cd $work_dir
# minimap2 original reads to TEs (?consensus reads)
for m in $(ls $input_dir/*.fq);do echo $m; filename=$(basename $m); minimap2 -ax map-ont  $TE_seqs $input_dir/$filename -t 2 --sam-hit-only --secondary=no -o 2TE_reads/2_${filename:0:-3}.sam;done

cd $work_dir/2TE_reads
## reads count for each chromsomes
for m in $(ls *.sam);do grep -v ^@ $m | cut -f1,3 | awk '{u=split($1,m,"_");print m[1], $2}' | sort | uniq| cut -d " " -f2 | sort | uniq -c | sort -nr > ${m:2:-4}.reads2chr.txt;done

## header modify and sort sam
for sam in $(ls 2_*.sam); do sed '/^@PG/d' $sam | sort -k1r| samtools view -S -b |samtools sort -@ 3 > 4_${sam:2:-4}.sorted.bam;done

## bam index
for m in $(ls 4_*.bam);do samtools index $m;done

#### write out TE-related reads
for m in $(ls $input_dir/*.fq);do echo $m; filename=$(basename $m); Rscript ../Fq.TE.extract.R ../$filename 4_${filename%.fq}.sorted.bam; done

## bam to bedgraph
for m in $(ls 4_*.bam);do bamCoverage -b $m --outFileFormat bedgraph -o 5_${m:2:-4}.bedgraph -p 3 --binSize 50; done

date



exit

## get sepcific region reads
samtools view  5_HCT116.modified.sorted.bam chr10:133779222-133780788 | cut -f1| awk '{u=split($1,m,"_");print m[1]}' | sort| uniq -c | sort >region_example/chr10:133779222-133780788.reads.txt
cat 2_chr10:133779222-133780788.reads.txt | while read LINE; do echo $LINE; cat /media/wuzefeng/文档/ZZ/200116.fq| awk 'BEGIN{RS="@";FS="\n"}; $1~/'$LINE'/{print ">"$1"\n"$2; exit}' >> 3_HCT_example_region.fa ;done

