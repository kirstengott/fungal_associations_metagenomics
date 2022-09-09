
import sys, os

fasta_in = sys.argv[1]
threads = sys.argv[2]
db_in = sys.argv[3]
output = sys.argv[4]

command = "diamond blastp --query {fasta} --db {db} --threads {p} --out {out_dir}".format(fasta = fasta_in, db = db_in, p = threads, out_dir = output)



if os.path.exists(output):
    print('ALREADY INITIATED:', command)
else:
    os.system(command)
    #print(command)
