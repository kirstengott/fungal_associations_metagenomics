
import os, sys, re, subprocess



sig_file = sys.argv[1]
db       = sys.argv[2]
out_dir  = sys.argv[3]

base = os.path.basename(sig_file)

bn = re.sub("\..*$", "", re.sub("_.*$", "", base))

db_base = os.path.splitext(os.path.basename(db))[0]

outfile = "{out}/{bn}_{db_base}.csv".format(out = out_dir, bn = bn, db_base = db_base)

command = "sourmash gather {sig} {db} --threshold-bp 10000 -o {out}".format(sig = sig_file, db = db, out = outfile)

checkpoint = outfile + "_chkpnt"

if os.path.exists(checkpoint):
    print('ALREADY DONE:', command)
else:
    #os.remove(checkpoint)
    print('Running:', command)
    p = subprocess.Popen(command, shell = True)
    os.waitpid(p.pid, 0)
    open(checkpoint, 'a').close()
