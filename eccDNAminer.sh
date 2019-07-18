mkidr pass
mkdir pass_merge

unzip pass.zip
for m in $(ls *.gz); do gunzip $m ; done
cat *.fastq > ../pass_merge/all.fastq
NanoPlot -t 40 --fastq  all.fastq --plots hex dot

## fastq convert fasta
sed -n '1~4s/^@/>/p;2~4p' fastq_runid_d9cccf88836673dacb0304fed8096fd4eca4f0fb_0.fastq > test.fasta

