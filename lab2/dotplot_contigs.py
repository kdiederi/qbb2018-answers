#!/usr/bin/env python3

"""
create dotplot of contig alignment to reference sequence (lastz data)
 
./dotplot_contigs.py <lastz_output.csv> <Sample/PlotName>

lastz_output from lastz generated at commandline using the output from
contig_compute.py for velvet or SPAdes
"""

import sys
import matplotlib.pyplot as plt


count = 0

for line in open(sys.argv[1]):
    fields = line.rstrip("\r\n").split()
    if line.startswith("#zstart1"):
        continue
    else:
        contig_length = int(fields[3]) - int(fields[2])
        plt.plot([int(fields[0]),int(fields[1])],[count,(count + contig_length)])
        count += contig_length

#reference (target) = zstart/end1, contigs (quary)= zstart/end2
#to set axis limits: plt.axis([xmin,xmax],[ymin,ymax])

plt.xlabel("Reference")
plt.xlim(0,110000)
plt.ylabel("Contigs")
plt.ylim(0,110000)
plt.title(str(sys.argv[2])+" Dotplot")
plt.savefig(str(sys.argv[2])+"_dotplot.png")
plt.close()