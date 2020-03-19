#!/bin/bash

#bash readme.sh
date

input_dir=/home/wuzefeng/Myresearch/ZZ/2020.3.19
work_dir=/home/wuzefeng/Myresearch/ZZ/2020.3.19
ref_genome=~/Myresearch/Genomes/DNA/GRCh38.primary_assembly.genome.fa
#faidx -x $ref_genome/GRCh38.primary_assembly.genome.fa # split chromsomes
ref_chromosomes=/home/wuzefeng/Myresearch/Genomes/DNA/human_chromsomes

# Tidehunter
for m in $(ls $input_dir/*.fq);do echo $m; filename=$(basename $m); TideHunter -o $work_dir/1_${filename%.fq}.cons.fa -t 3 $input_dir/$filename;done

cd work_dir
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

## bam to bedgraph
for m in $(ls 5_*.bam);do bamCoverage -b $m --outFileFormat bedgraph -o 6_${m:2:-4}.bedgraph -p 3 --binSize 50; done

date

