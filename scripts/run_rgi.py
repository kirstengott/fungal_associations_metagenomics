


import sys,os

fasta = sys.argv[1]

output_dir = os.path.dirname(fasta)
outfile    = "{}/card_out".format(output_dir)




fasta_temp          = fasta + ".card_temp.faa"
clean_fasta_command =  'cat {} | sed -e "s/*//" > {}'.format(fasta, fasta_temp)
os.system(clean_fasta_command)


cmd        = "rgi main -i {faa} -o {out} -a DIAMOND -n 1 --include_loose --clean -t protein".format(faa = fasta_temp, out = outfile)








if os.path.exists(outfile + ".txt"):
    print('ALREADY INITIATED:', cmd)
else:
    #print(cmd)
    print(fasta)
    os.system(cmd)


os.remove(fasta_temp)
