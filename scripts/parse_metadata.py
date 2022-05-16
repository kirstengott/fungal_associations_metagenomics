import os
import re, sys
from config  import *


keys_keep = ['JGI_DIR_ID',
             'pipeline_version',
             'Proposal_Name',
             'Proposal_ID',
             'Analysis_Project/Task_ID',
             'Sequencing_Project_ID(s)',
             'Library',
             'Platform',
             'FilteredData',
             'RawData',
             'The_number_of_raw_input_reads_is',
             'The_number_of_reads_filtered_out_due_to_artifact_or_low_quality_is',
             'A',
             'C',
             'G',
             'T',
             'N',
             'IUPAC',
             'Other',
             'GC',
             'GC_stdev',
             'Main_genome_scaffold_total',
             'Main_genome_contig_total',
             'Main_genome_scaffold_sequence_total',
             'Main_genome_contig_sequence_total',
             'Main_genome_scaffold_N/L50',
             'Main_genome_contig_N/L50',
             'Main_genome_scaffold_N/L90',
             'Main_genome_contig_N/L90',
             'Max_scaffold_length',
             'Max_contig_length',
             'Number_of_scaffolds_>_50_KB',
             '%_main_genome_in_scaffolds_>_50_KB',
             'The_number_of_reads_used_as_input_to_aligner_is',
             'The_number_of_aligned_reads_is',
             'The_number_of_unmapped_read_pairs_as_input_to_bbmerge_is',
             'The_number_of_merged_read_pairs_is',
             'If_you_have_any_questions,_please_let_us_know',
             'GOLD_ID_(BIOSAMPLE)',
             'BIOME',
             'FEATURE',	
             'MATERIAL',
             'LATITUDE_AND_LONGITUDE',
             'VERTICAL_DISTANCE_ALTITUDE',
             'GEOGRAPHIC_LOCATION',
             'SAMPLE_COLLECTION_DATE',
             'IMG_ID_(IMG_TAXON_OID)',
             'STUDY',
             'GOLD_ID_(SEQUENCING_PROJECT)',
             'GOLD_ID_(ANALYSIS_PROJECT)',
             'NCBI_BIOPROJECT']

def parse_pipeline_version(file_path):
    with open(file_path, 'r') as fh:
        version = ''
        for line in fh:
            items = line.strip().split(" ")
            version = items[2]
            break
        return(version)

def parse_readme(file_path):
    iupac_flag = 0
    with open(file_path, 'r') as fh:
        out_dict = {}
        for line in fh:
            if iupac_flag == 1:
                values = line.strip().split()
                assemb_stats = dict(zip(keys, values))
                out_dict.update(assemb_stats)
                iupac_flag = 0
            if ':' in line:
                items = line.strip().split(":")
                items = [x.strip() for x in items]
                items = [x for x in items if x]
                items = [re.sub(" ", "_", x) for x in items]
                if len(items) == 2:
                    out_dict[items[0]] = items[1]
            if 'IUPAC' in line:
                keys = line.strip().split()
                iupac_flag = 1
        return(out_dict)


def parse_study_info(file_path):
    dict_out = {}
    with open(file_path, 'r') as fh:
        i = 0
        for line in fh:
            if i == 0:
                i = 1
                pass
            else:
                items = re.sub(" ", "_", line.strip())
                items = re.split('__+', items)
                if len(items) == 2:
                    val = items[1]
                    key = items[0]
                    dict_out[key] = val
                else:
                    if len(items) == 1:
                        continue
                    else:
                        print('DICTIONARY CONSTRUCTION ERROR AT "def parse_study_info"', file_path, items)
                    
    return(dict_out)
                
            
    
file_patterns = ['README.txt', 'pipeline_version.info', 'Study_Information', 'Sample_Information']

metadata_table_fh = open(metadata_table, 'w')

key_init = 0
print_line = 0

for i in os.listdir(data_dir):
    dir_path = data_dir + i
    if not os.path.isdir(dir_path):
        continue
    dir_files = os.listdir(dir_path)
    data_in = {}
    ## iterate through the files in each directory
    for file_i in dir_files:
        data_in['JGI_DIR_ID'] = i
        ## iterated through the patterns I want to match
        for file_p in file_patterns:
            ## if  I match a file that I want, parse the data
            if file_p in file_i:
                file_path_in = dir_path + "/" + file_i

                ### parse the first file
                if file_p == file_patterns[0]:
                    readme_dict = parse_readme(file_path_in)
                    data_in.update(readme_dict)

                ## parse the second file
                if file_p == file_patterns[1]:
                    pipeline_version = parse_pipeline_version(file_path_in)
                    data_in['pipeline_version'] = pipeline_version
                ## parse the third file
                if file_p == file_patterns[2]:
                    study_info = parse_study_info(file_path_in)
                    # study_info.pop('SRA')
                    # study_info.pop('')
                    data_in.update(study_info)
                ## parse the fourth file
                if file_p == file_patterns[3]:
                    sample_info = parse_study_info(file_path_in)
                    data_in.update(sample_info)
                
            else:
                continue
    # if print_line == 100:
    #     for d in data_in:
    #         print(d, "\t", data_in[d])
    #     sys.exit()
    # else:
    #     print_line += 1
    
    if key_init == 0:
        key_init = 1
        metadata_table_fh.write("\t".join(keys_keep) + "\n")
    values_out = []
    for i in keys_keep:
        if i in data_in:
            values_out.append(data_in[i])
        else:
            values_out.append('NA')
    metadata_table_fh.write("\t".join(values_out) + "\n")
    data_in = {}

metadata_table_fh.close()
