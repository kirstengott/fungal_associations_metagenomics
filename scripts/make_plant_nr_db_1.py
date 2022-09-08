import os, sys, re

def listToDict(lst):
    op = { i : 1 for i in lst }
    return op

## first parse in the plant accessions

plant_taxa_ids = 'db/blast/taxonomy_plant.txt'

plant_taxa = []

with open(plant_taxa_ids, 'r') as fh:
    # plant_nr = [re.sub("\..*$", "", line.strip()) for line in fh]
    plant_nr = [line.strip() for line in fh]


## next parse in the matching taxids





plant_taxa_d = listToDict(plant_taxa)


prot2taxid = 'db/blast/prot.accession2taxid.FULL' ## column 1 is the accession, column 2 is the taxid


prot2tax = {}

print('First Done')

count = 0
with open(prot2taxid, 'r') as fh:
    for line in fh:
        if count == 0:
            count = 1
            continue
        else:
            items = line.strip().split()
            if items[1] in plant_taxa_d.keys():
                prot2tax[items[0]] = items[1]
            else:
                continue


print('Second Done')
print(prot2taxid[prot2taxid.keys()[0]])


## filter the nr fasta and rename them according to the kaiju spec (counter_taxid)

nr       = 'db/blast/nr.fa'


file_out = 'db/blast/nr_plant_2.fa'
fo       = open(file_out, 'w')


counter  = 1
seq_flag = 0
with  open(nr, 'r') as fh:
    for line in fh:
        if seq_flag == 1:
            fo.write(line)
        if line.startswith('>'):
            seq_flag = 0
            acc = re.sub(">", "", line.strip().split()[0])
            if acc in prot2tax.keys():
                new_fa_id = ">" + str(counter) + "_" + prot2tax[acc] + "\n"
                fo.write(new_fa_id)
                counter += 1
                seq_flag = 1
        else:
            continue
    





# Custom database
# It is also possible to make a custom database from a collection of protein sequences. The format needs to be a FASTA file in which the headers are the numeric NCBI taxon identifiers of the protein sequences, which can optionally be prefixed by another identifier (e.g. a counter) followed by an underscore, for example:

# >1_1358
# MAQQRRGGFKRRKKVDFIAANKIEVVDYKDTELLKRFISERGKILPRRVTGTSAKNQRKVVNAIKRARVMALLPFVAEDQN
# >2_44689
# MASTQNIVEEVQKMLDTYDTNKDGEITKAEAVEYFKGKKAFNPERSAIYLFQVYDKDNDGKITIKELAGDIDFDKALKEYKEKQAKSKQQEAEVEEDIEAFILRHNKDDNTDITKDELIQGFKETGAKDPEKSANFILTEMDTNKDGTITVKELRVYYQKVQKLLNPDQ
# >3_352472
# MKTKSSNNIKKIYYISSILVGIYLCWQIIIQIIFLMDNSIAILEAIGMVVFISVYSLAVAINGWILVGRMKKSSKKAQYEDFYKKMILKSKILLSTIIIVIIVVVVQDIVINFILPQNPQPYVYMIISNFIVGIADSFQMIMVIFVMGELSFKNYFKFKRIEKQKNHIVIGGSSLNSLPVSLPTVKSNESNESNTISINSENNNSKVSTDDTINNVM
# >4_91061
# MTNPFENDNYTYKVLKNEEGQYSLWPAFLDVPIGWNVVHKEASRNDCLQYVENNWEDLNPKSNQVGKKILVGKR

