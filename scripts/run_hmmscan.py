



import os,sys

fasta = sys.argv[1]
db = sys.argv[2]
threads = sys.argv[3]


output_dir = os.path.dirname(fasta)
outfile    = "{}/pfam.out".format(output_dir)



cmd = "hmmscan -o {out} --tblout {output_dir}/seq_pfam.table --domtblout {output_dir}/dom_pfam.table -E .01 --domE .01 --cpu {p} {db} {input}".format(out = outfile, output_dir = output_dir, p = threads, db = db, input = fasta)



if os.path.exists(outfile):
    print('ALREADY INITIATED:', cmd)
else:
    os.system(cmd)
#    print(cmd)
