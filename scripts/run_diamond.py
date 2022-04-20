
import sys, os

fasta_in = sys.argv[1]
threads = sys.argv[2]
db_in = sys.argv[3]
output_dir = os.path.dirname(fasta_in)

command = "diamond blastp --query {fasta} --db {db} --threads {p} --out {out_dir}/diamond_blastp.tsv".format(fasta = fasta_in, db = db_in, p = threads, out_dir = output_dir)


outfile = "{out_dir}/diamond_blastp.tsv".format(out_dir = output_dir)


if os.path.exists(outfile):
    print('ALREADY INITIATED:', command)
else:
    #os.system(command)
    print(command)
