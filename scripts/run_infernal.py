
import sys
import os
import subprocess

fasta = sys.argv[1]



output_dir = os.path.dirname(fasta)
prefix = os.path.dirname(fasta).split("/")[-1] + os.path.basename(fasta)


bp_com = "esl-seqstat {}"
output = subprocess.check_output(['esl-seqstat', fasta]).decode(sys.stdout.encoding)
out = output.split("\n")

total_n = ''

for x in out:
    if 'Total' in x:
        total_n = x.split(' ' )
    else:
        pass

total_bp = int(total_n[-1])/1000000



infernal_cmd = "cmscan -Z {bp} --cpu 1 --cut_ga --rfam --nohmmonly --tblout {out_dir}/{pre}_infernal-genome.tblout --fmt 2 --clanin ~/db/Rfam.clanin ~/db/Rfam.cm {fasta_in}  > {out_dir}/{pre}_infernal-genome.cmscan".format(bp = total_bp, pre = prefix, fasta_in = fasta, out_dir = output_dir)





