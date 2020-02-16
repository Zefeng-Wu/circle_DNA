#!/bin/python
import sys
infile = sys.argv[1]
outfile = sys.argv[2]
outfile = open(outfile,"w")
with open(infile) as fh:
    for row in fh:
        if row.startswith("@"):
            outfile.write(row)
        else:
            data = row.strip().split("\t")
            #print (data)
            data[0] = data[0].split("_")[0]+"_"+data[0].split("_")[1]
            outfile.write("\t".join(data)+"\n")
print ("OK!")
