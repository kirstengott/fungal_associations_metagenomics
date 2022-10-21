
import os, sys, re, subprocess



sig_file = sys.argv[1]
#db       = sys.argv[2]
out_dir  = sys.argv[2]

base = os.path.basename(sig_file)

bn = re.sub("\..*$", "", re.sub("_.*$", "", base))

#db_base = os.path.splitext(os.path.basename(db))[0]
db = "/home/gotting/projects/fungal_associations_metagenomics/db/sourmash/genbank-2022.03-archaea-k31.zip /home/gotting/projects/fungal_associations_metagenomics/db/sourmash/genbank-2022.03-bacteria-k31.sbt.zip /home/gotting/projects/fungal_associations_metagenomics/db/sourmash/genbank-2022.03-fungi-k31.zip /home/gotting/projects/fungal_associations_metagenomics/db/sourmash/genbank-2022.03-protozoa-k31.zip /home/gotting/projects/fungal_associations_metagenomics/db/sourmash/genbank-2022.03-viral-k31.zip /home/gotting/projects/fungal_associations_metagenomics/db/sourmash/genbank-2022.08-plant-k31.sbt.zip"

outfile = "{out}/{bn}.csv".format(out = out_dir, bn = bn)

command = "sourmash gather {sig} {db} --threshold-bp 10000 -o {out} 2> {out}_error.txt 1> {out}_screen_output.txt".format(sig = sig_file, db = db, out = outfile)

checkpoint = outfile + "_chkpnt"

if os.path.exists(checkpoint) and os.path.exists(outfile):
    print('ALREADY DONE:', command)
else:
    #os.remove(checkpoint)
    print('Running:', command)
    p = subprocess.Popen(command, shell = True)
    os.waitpid(p.pid, 0)
    open(checkpoint, 'a').close()



