


import sys, os, re


sam = sys.argv[1]
threads = sys.argv[2]

outfile = re.sub('.sam.gz', '.bam', sam)
sam_f = re.sub('.gz', '', sam)

cmd = "gunzip {samfile}.gz; samtools sort -@ {p} -m 10G {samfile} -o {outfile}; gzip {samfile}".format(p = threads, samfile = sam_f, outfile = outfile)


if os.path.exists(outfile):
    print('ALREADY INITIATED:', cmd)
else:
    os.system(cmd)
#    print(cmd)


