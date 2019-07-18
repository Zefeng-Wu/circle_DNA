mkdir pass_merge
unzip pass.zip
for m in $(ls *.gz); do gunzip $m ; done
cat *.fastq > ../pass_merge/all.fastq

### reads quality statistics
NanoPlot -t 40 --fastq  all.fastq --plots hex dot

## fastq convert fasta
sed -n '1~4s/^@/>/p;2~4p' all.fastq > all.fasta

## perform tandem repeat finder to seach TRs in reads
~/soft/trf409.linux64 all.fasta 2 7 7 80 10 50 2000 -l 6 -h > out.txt
# match score: 2
# mismatch score: 7
# indel score: 7
# acth probality: 80 
# min alignemt score to report: 50
# max period: 2000
# -l longest repeat size(million) 

### format convert
python test.py in.dat out.txt

