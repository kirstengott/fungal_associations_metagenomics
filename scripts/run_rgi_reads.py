### TODO

## 1) Add code for paired end reads
## 2) Add code to gunzip -c multiple fastqs and pipe into the command


import sys,os
from pathlib import Path




fastq   = sys.argv[1] ## accepted gzipped fastqs
out_pre = sys.argv[2]
threads = sys.argv[3]



touch = out_pre + '.chkpnt'




if len(sys.argv) > 4 :
    cmd        = "rgi bwt -1 {fastq} -2 {fq2} -a bwa -n {p} -o {outfile} --clean ".format(fastq = fastq, outfile = out_pre, p = threads, fq2 = sys.argv[4])
else:
    cmd        = "rgi bwt -1 {fastq} -a bwa -n {p} -o {outfile} --clean ".format(fastq = fastq, outfile = out_pre, p = threads)



if os.path.exists(touch):
    print('ALREADY INITIATED:', cmd)
else:
    Path(touch).touch()
    os.system(cmd)

