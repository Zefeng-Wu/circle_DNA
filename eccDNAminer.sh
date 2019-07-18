mkidr pass
mkdir pass_merge

unzip pass.zip
for m in $(ls *.gz); do gunzip $m ; done
cat *.fastq > ../pass_merge/all.fastq
NanoPlot -t 40 --fastq  all.fastq --plots hex dot

## fastq convert fasta
sed -n '1~4s/^@/>/p;2~4p' fastq_runid_d9cccf88836673dacb0304fed8096fd4eca4f0fb_0.fastq > test.fasta



./trf409.linux64  /home/wuzefeng/MyResearch/genome.db/TAIR/dna/Arabidopsis_thaliana.TAIR10.31.dna.toplevel.fa 2 7 7 80 10 100 2000 -d -h

./trf409.linux64 fasta文件名 match得分 错配得分 indel罚分 PM PI 最小比对得分 tandem覆盖的序列长度 -d -h(阻止html)


