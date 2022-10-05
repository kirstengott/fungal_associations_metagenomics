import sys, os


rf     = ''
scaf   = ''
start  = ''
end    = ''
strand = ''
        
out = "{},{},{},{},{},{}"

out_flag = 0


infile = sys.argv[1]

sort_cmnd = "cat {inf} | tr -s ' ' | tr ' ' , | sort -t, -k4,4n -k10,10n -k11,11n >>{inf}_sorted.csv".format(inf = infile)

os.system(sort_cmnd)

tmp_file = '{}_sorted.csv'.format(infile)

with open(tmp_file, 'r') as fh:
    for line in fh:
        if '#' in line:
            continue
        else:
            out_flag = 0
            liner = line.strip().split(',')
            if len(liner) < 17:
                continue
            if rf == '':
                rf     = liner[2]
                scaf   = liner[3]
                start  = int(liner[9])
                end    = int(liner[10])
                strand = liner[11]
                evalue  =  float(liner[17])
                continue
            else:
                rf_n     = liner[2]
                scaf_n   = liner[3]
                start_n  = int(liner[9])
                end_n    = int(liner[10])
                strand_n = liner[11]
                evalue_n  =  float(liner[17])
            if scaf != scaf_n:
                out_flag = 1
            else:
                if strand_n == strand: ## check strand
                    if strand == "+":
                        start_h = start - 300 ## set the bounds to be the length of a small gene
                        end_h   = end + 300
                        if start_n >= start_h and end_n <= end_h: ## keep the one with the smaller eval
                            if evalue_n < evalue:
                                rf     = rf_n
                                scaf   = scaf_n
                                start  = start_n
                                end    = end_n
                                strand = strand_n
                                evalue = evalue_n
                            continue
                        else: ## features next to eachother, keep both
                            out_flag = 1
                    else:     ## logic for '-' strand
                        start_h = start + 300 ## set the bounds to be the length of a small gene
                        end_h   = end - 300
                        if start_n <= start_h and end_n >= end_h: ## only keep the first
                            if evalue_n < evalue:
                                rf     = rf_n
                                scaf   = scaf_n
                                start  = start_n
                                end    = end_n
                                strand = strand_n
                                evalue = evalue_n
                            continue
                        else: ## features next to eachother, keep both
                            out_flag = 1
                else:
                    out_flag = 1 ## features on opposite strands
        if out_flag == 1:
            print(out.format(rf, scaf, start, end, strand, evalue)) ## print the first and reset to next
            rf     = rf_n
            scaf   = scaf_n
            start  = start_n
            end    = end_n
            strand = strand_n
            evalue = evalue_n
            
print(out.format(rf, scaf, start, end, strand, evalue)) ## print the first and reset to next

                

                
