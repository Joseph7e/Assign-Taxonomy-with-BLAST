#!/mnt/lustre/software/linuxbrew/colsa/bin/python3
#Author: Joseph Sevigny
#Affiliation: Hubbard Center for Genome Studies, University of New Hampshire
#Date: 03/15/2017

import sys, re, os, subprocess, argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#OPTIONAL ARGUMENTS
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
#parser.add_argument("-h","--help","-help", help="prints help and exits")
# Consensus capture options
parser.add_argument("--cutoff_species", help="cutoff for finest taxonomic level", type=int, default=97)
parser.add_argument("--cutoff_family", help="cutoff for family taxonomic level", type=int, default=90)
parser.add_argument("--cutoff_phylum", help="cutoff for phylum taxonomic level, also acts as ultimate cutoff value for blast", type=int, default=80)

# Filtering options
parser.add_argument("--length_percentage", help="cutoff for query_hit/length_of_query i.e. query coverage", type=float, default=0.0)
parser.add_argument("--length_cutoff", help="primary cutoff for length of hit", type=int, default=0)
parser.add_argument("--hits_to_consider", help="number of hits to consider when gathering consensus taxonomy", type=float, default=3)
parser.add_argument("--percent_sway", help="when comparing greater than 1 blast hit, value is the percent from best hit considered. i.e. if best blast hit is 99.5 ID, a value of 0.5 will consider everything 99 and greater when creating the consensus", type=float, default=0.5)

# BLAST options
parser.add_argument("--blast_evalue", help="setting for e-value cutoff for blast, must be in form 1e-X", type=str, default="1e-10")
parser.add_argument("--blast_threads", help="set the number of threads for blast", type=str, default="24")
parser.add_argument("--blast_flavor", help="select blastp, tblastn etc, must have correct database and query sequence formats", type=str, default="blastn")
parser.add_argument("--blast_database", help="path tp blast database, be sure to run makeblastdb on the database first, type 'IGNORE' if precomputed blast is given")
parser.add_argument("--blast_file", help="precomputed blast results, MUST BE MY CUSTOMIZED FORMAT, '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore'")

# Other options
parser.add_argument("--ncbi_nt", help="REQUIRED flag for use of ncbi nt database", action="store_true")
parser.add_argument("--output_dir", help="output directory name", type=str, default="Assigned_Taxonomy")
parser.add_argument("--config_file", help="NOT YET IMPLEMENTED, an optional parameters file, any arguments in this file will override other options")

#REQUIRED ARGUMENTS
parser.add_argument("sequence_file", help="seqs.fna file from qiime or any multifasta, just make sure header has unique id with a space")

parser.add_argument("tax_file", help="path to silva or customized blast database taxonomy file")

#parser.add_argument
args = parser.parse_args()

if args.blast_file == None and args.blast_database == None:
    print ('The argument --blast_file or --blast_database is required. You must provide at least one.')
    sys.exit()

##### CONSTRUCT OUTTDIR NAME
temp_outdir = args.output_dir
if temp_outdir[-1] == ('/'):
    temp_outdir = temp_outdir[:-1]
current_directory_files = os.listdir('./')

at_count = 0
outdir = temp_outdir
while True:
    if outdir in current_directory_files:
        at_count += 1
        outdir = temp_outdir+'_'+str(at_count)
    else:
        outdir += '/'
        break

#OUTPUT_FILES
os.mkdir(outdir)
taxonomy_assignment_outfile = outdir + 'taxonomy_assignment_per_sequence.tsv'
log_file_handle = open(outdir + 'log_file.txt','w')
taxonomy_assignment_seven_levels_outfile = outdir + 'taxonomy_assignment_per_sequence_eight_levels.tsv'
#config_out_handle = open(outdir + 'config_file.txt','w')

## DONE

def log_and_print(statement, log_file_name_handle=log_file_handle, print_bool=args.verbose):
    ''' takes care of log creation and print statements'''
    if print_bool == True:
        print (statement)
    log_file_name_handle.writelines(statement+'\n')





log_and_print('####VERBOSE OPTION GIVEN\n')
log_and_print('#MAIN ARGUMENTS')
log_and_print('sequence_file: ' + args.sequence_file)
log_and_print('blast_database: ' + str(args.blast_database))
log_and_print('taxonomy_file: '+ args.tax_file)
log_and_print('\n#Filtering OPTIONS')
log_and_print('blast_hits_to_consider: ' + str(args.hits_to_consider))
log_and_print('sway_from_best_hit: ' +  str(args.percent_sway))
log_and_print('phylum_level_cutoff: ' +  str(args.cutoff_phylum))
log_and_print('family_level_cutoff: ' +  str(args.cutoff_family))
log_and_print('species_level_cutoff: ' +  str(args.cutoff_species))
log_and_print('length_of_hit_cutoff: ' +  str(args.length_cutoff))
log_and_print('proportion_of_query_coverage: ' +  str(args.length_percentage))
log_and_print('blast e-value cutoff: ' +  args.blast_evalue)
log_and_print('\nOther OPTIONS')
log_and_print('BLAST_file given?: ' +  str(args.blast_file))
log_and_print('BLAST flavor: ' +  args.blast_flavor)
log_and_print('BLAST evalue cutoff: ' +  args.blast_evalue)
log_and_print('BLAST threads: ' +  args.blast_threads)
log_and_print('NCBI flag?: ' +  str(args.ncbi_nt))
log_and_print('Output dir: ' +  outdir)
log_and_print('Param file: ' +  str(args.config_file))



log_and_print('\n###Program Start\n')


###INFO FOR SILVA 128 99 OTUS all levels
taxonomy_categories = ['superkingdom', 'subkingdom', 'sub_subkingdom', 'kingdom', 'tmp1', 'tmp2', 'phylum', 'class', 'family', 'genus', 'species', 'tmp3', 'tmp4','tmp5', 'tmp6']

#ITS database
#taxonomy_categories = ['kingdom','phylum','class','order','family','genus','species']


# 16S database
### OTHER OPTIONS BASED ON DATABASE OF CHOICE
uncultured_cutoff_level = len(taxonomy_categories)
species_level = 14
family_level = 10
phylum_level = 7
# #taxonomy_to_grab = [0,3,6,10,11,12] # file for qiime output
taxonomy_to_grab = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]


#### ITS DATABASE
#species_level = 6
#family_level = 4
#phylum_level = 1

#taxonomy_to_grab = [0,1,2,3,4,5,6]



if args.ncbi_nt:
    ###INFO FOR TAXONOMY DATABASES
    taxonomy_categories = ['superkingdom', 'kingdom', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'superorder', 'order', 'superfamily', 'family', 'subfamily', 'genus', 'species']
    #qiime_all_level_12_categories = [superkingdom, subkingdom, sub_subkingdom, kingdom, tmp1, tmp2, phylum, subphylum, class, order, family, genus, species]
    seven_levels = [0,1,2,5,8,10,12,13]
    ### OTHER OPTIONS BASED ON DATABASE OF CHOICE
    uncultured_cutoff_level = len(taxonomy_categories)
    species_level = 13
    family_level = 10
    phylum_level = 2

    taxonomy_to_grab = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]

else:
    seven_levels = [0,1,2,5,8,10,12,13]



log_and_print("#Taxonomy Categories\n"+':'.join(taxonomy_categories))
log_and_print("#The below levels should match the target, i.e. species_level = species")
log_and_print("Species_level =  " + str(species_level) + ' = ' + taxonomy_categories[species_level])
log_and_print("Family_level =  " + str(family_level) + ' = ' + taxonomy_categories[family_level])
log_and_print("Phylum_level =  " + str(phylum_level) + ' = ' + taxonomy_categories[phylum_level] +'\n')




#Potential Additions:
## --> Add option to blast more than one otu if they might be different, make the otu selection strigent



#####Dictionaries

sequence_dict = SeqIO.to_dict(SeqIO.parse(args.sequence_file, "fasta"))
taxonomy_dictionary = {} # tax_code:[taxonomy]
sequence_taxonomy_dict = {} # header:[best_taxonomy]
percent_id_dict = {}






###do_blast
def do_blast(query, database, blast_output):

    '''performs blasts and returns the results'''
    #custom_blast_format = '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore slen staxids'
    custom_blast_format = '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore slen staxids'

    #db_command = ["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", "temp_db"] #command to 1econstruct database
    #sd = subprocess.Popen(db_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE); sd.communicate()
    blast_command = [args.blast_flavor, "-query", query, "-db", database , "-num_threads", args.blast_threads, "-evalue", args.blast_evalue, "-out", blast_output, "-outfmt", custom_blast_format] #command to complete blastn
    log_and_print('#BLAST COMMAND')
    log_and_print(' '.join(blast_command))
    sp = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)#blast = sp.communicate()
    err = sp.stderr.read()
    if err:
        log_and_print(err)
        log_and_print('ERROR: with blast above')
        sys.exit()
    sp.communicate()

if args.blast_file == None:
    log_and_print('Running blast on input sequences (this step takes the longest, maybe run overnight)')
    blast_file = outdir + 'blast_output_custom_format.tsv'
    do_blast(args.sequence_file,args.blast_database,blast_file)
else:
    blast_file = args.blast_file



#parse taxonomy file

log_and_print('Parsing taxonomy file and storing info in dictionary (might take a moment if using nt database)')
for line in open(args.tax_file,'r'):
    elements = line.rstrip().split('\t')
    tax_id = elements[0]
    tax_info =re.sub('D_[0-9]*__','',elements[1]).replace(' ','_').split(';') # reformat qiime tax file
    taxonomy_dictionary[tax_id]=tax_info
    #taxonomy_dictionary[tax_id]=[tax_info[3],tax_info[6:]] ???



####PARSE BLAST and ASSIGN TAXONOMY#######
#make so assignment avoids one mis when the rest agree


def Assign_Taxonomy(current_query, current_best_hits):
    ''' the leg work'''
    log_and_print('\nASSIGNING TAXONOMY FOR ' + current_query + ' total hits passing initial filters = ' + str(len(current_best_hits)))
    ##FILTER FOR THE BEST HIT and set cutoff

    #find best hit based on percent id
    top_hits = [] #should add check for best bitscore
    for key, value in current_best_hits.items():#TURN INTO LAMBDA
        top_hits.append(float(value[0]))

    #Fill best_blast_dictionary with those within sway distance of best blast hit
    best_blast_dictionary = {} #unique_blast_id:[percent_id, tax_code]
    for key, value in current_best_hits.items():
        if float(value[0]) > (max(top_hits)) - args.percent_sway : # *****important cutoff
            best_blast_dictionary[key]=value
            log_and_print(value[1] + ' ' + value[0] +  ' --> CAPTURED after percent sway filter')
    ##Assign Taxonomic Certainty
    level = phylum_level #'phylum'
    if max(top_hits) >= args.cutoff_family:
        level = family_level #'family'
    if max(top_hits) >= args.cutoff_species:
        level = species_level #'species'

    log_and_print("## Providing consensus taxonomy up to level " + str(level) + " : " + taxonomy_categories[level] )

    blast_percent = max(top_hits)

    ##Assign Taxonomy
    temp_taxonomy_dict = {} #unique_code:[taxonomy]
    count = 0
    for key,value in best_blast_dictionary.items():
        count += 1
        pkey = str(value[0])+ '_'+str(count)
        #Mask uncultured taxonomy
        if value[1] not in taxonomy_dictionary.keys(): # make sure that the key is in the dictionary
            taxonomy_dictionary[value[1]] = taxonomy_categories
            pkey = 'X'+pkey
        if 'uncultured' in ' '.join(taxonomy_dictionary[value[1]][:uncultured_cutoff_level+1]) or 'unidentified' in ' '.join(taxonomy_dictionary[value[1]][:uncultured_cutoff_level+1]):
            pkey = 'X'+pkey
            #print ('UNCULTURED', current_query)
        temp_taxonomy_dict[pkey]=taxonomy_dictionary[value[1]]#[:level]
        log_and_print(pkey + ' ' + ';'.join(taxonomy_dictionary[value[1]]))


    #Determine best common taxonomy up to set certainty (certainty = cutoff level)
    # for l, t in temp_taxonomy_dict.items():
    #     print (l, t)

    best_level_taxonomy = [] #final taxonomy assignment
    set_list = []
    for i in range(level+1): # start at kingdom and move toward species until difference is met.
        s = set() # set to hold labels from all tax code.
        the_okay = False

        #avoid error when all sequences are uncultured
        for one in temp_taxonomy_dict.keys():
            if not one.startswith('X'):
                the_okay = True

        for k,j in temp_taxonomy_dict.items():
            if not k.startswith('X') and len(temp_taxonomy_dict) > 1 and the_okay == True:
                s.add(j[i])
            if the_okay == False:
                s.add(j[i])
        set_list.append(s)
        if len(s) == 1:
            e = next(iter(s))
            best_level_taxonomy.append(e)
    #print (set_list)
    #when there is only one hit or they all match the same thing
    if not bool(best_level_taxonomy):
        for value in temp_taxonomy_dict.values():
            best_level_taxonomy = value[:level+1]

    best_level_taxonomy = list(filter(bool, best_level_taxonomy))#remove empty taxonomic levels
    needed_unknowns = (len(taxonomy_categories)-len(best_level_taxonomy)) # Add undetermined
    for i in range(needed_unknowns):
        best_level_taxonomy.append('undetermined')

    ### FINALLY FILL OTU INFORMATION!!!!
    return best_level_taxonomy, blast_percent



####################################
log_and_print('Parsing Blast and Assigning Taxonomy\n\n########################')

current_best_hits = {} #count_code:[percent_id, tax_code]
current_query = ''

#total counts
total_sequences_in_blast = 0
total_sequences_assigned = 0


top_hit_dict = {}

for line in open(blast_file,'r'):
    elements = line.rstrip().split('\t'); sequence_id = elements[0]; tax_code = elements[2]
    percent_id = elements[3]; length_hit = elements[4]; length_query = elements[1]; ncbi_id = elements[-1]
    if sequence_id not in top_hit_dict.keys():
        top_hit_dict[sequence_id] = [percent_id,length_hit,length_query]
    if args.ncbi_nt:
        if ';' in ncbi_id:
            ncbi_id = ncbi_id.split(';')[0]
        tax_code = ncbi_id

    if current_query != sequence_id: #INITIATION of new query

        #turn below into a function or definition
        total_sequences_in_blast += 1
        if bool(current_best_hits) and current_query: #ENSURE A SET OF BLASTS TO PARSE

            best_level_taxonomy, blast_percent = Assign_Taxonomy(current_query, current_best_hits)

            # ### FINALLY FILL Sequence INFORMATION!!!!
            total_sequences_assigned += 1
            sequence_taxonomy_dict[current_query] = best_level_taxonomy
            percent_id_dict[current_query] = blast_percent
            log_and_print('\nTaxonomy Assignment for ' +  current_query +  ' = ' +  ':'.join(best_level_taxonomy)+'\n\n\n######')

        #reset query_id
        current_query = sequence_id
        current_best_hits = {}

    ## Filter Hits and Loop Through current Query and add up to x number of blast hits
    if int(length_hit)/int(length_query)>args.length_percentage and float(percent_id)>args.cutoff_phylum and int(length_hit)>args.length_cutoff:
        if len(current_best_hits.keys()) < args.hits_to_consider:
            log_and_print('#BLAST LINE : ' + line.rstrip())
            label = len(current_best_hits.keys()) # provide unique key based on number of matched
            current_best_hits[label] = [percent_id, tax_code]
####################################################
### Get the last hits info
best_level_taxonomy, blast_percent = Assign_Taxonomy(current_query, current_best_hits)

# ### FINALLY FILL Sequence INFORMATION!!!!
total_sequences_assigned += 1
sequence_taxonomy_dict[current_query] = best_level_taxonomy
percent_id_dict[current_query] = blast_percent
log_and_print('Taxonomy Assignment for ' +  current_query +  ' = ' +  ':'.join(best_level_taxonomy)+'\n\n\n######')


### QIIME FORMATTED OUTPUT

if args.verbose:
    print ('\nConstructing output file')

#Write output

for seq in sequence_dict.keys():
    seven_levels_out_line = seq + ' '
    output_line = seq + ' '
    if seq in sequence_taxonomy_dict.keys():
        level_count = 0
        for taxonomy_level in sequence_taxonomy_dict[seq]:
            if level_count in seven_levels:
                seven_levels_out_line += (taxonomy_level + ';')
            output_line += (taxonomy_level + ';')
            level_count += 1
        output_line = output_line.rstrip(';')
        seven_levels_out_line = seven_levels_out_line.rstrip(';')
        with open(taxonomy_assignment_outfile,'a') as t:
            #t.writelines(output_line+ ' ' +str(percent_id_dict[seq]) +'\n')
            t.writelines(output_line+ ' ' +':'.join(top_hit_dict[seq]) +'\n')
        with open(taxonomy_assignment_seven_levels_outfile,'a') as t:
            #t.writelines(seven_levels_out_line+ ' ' +str(percent_id_dict[seq]) +'\n')
            t.writelines(seven_levels_out_line+ ' ' +':'.join(top_hit_dict[seq]) +'\n')
    else: # handle cases that were not in blast file
        output_line += ('Unassigned;'*(len(taxonomy_to_grab)-1)) + 'Unassigned 0'
        seven_levels_out_line += ('Unassigned;'*7) + 'Unassigned 0'
        with open(taxonomy_assignment_outfile,'a') as t:
            t.writelines(output_line+'\n')
        with open(taxonomy_assignment_seven_levels_outfile,'a') as t:
            t.writelines(seven_levels_out_line+'\n')

total_sequences = len(sequence_dict.keys())

log_and_print('Total number of sequences in fasta: ' +  str(total_sequences))
log_and_print('Total number of sequences in blast_file: ' +  str(total_sequences_in_blast))
log_and_print('Total number sequences assigned to taxonomy: ' +  str(total_sequences_assigned))
log_and_print('Percentage sequences with assigned taxonomy: ' + str(total_sequences_assigned/total_sequences * 100)   + '%')


## sort the taxonomy file for viewing
os.system("sort -k2 " + outdir + "taxonomy_assignment_per_sequence.tsv > " + outdir + "tmp.txt && mv " + outdir + "tmp.txt " + outdir + "taxonomy_assignment_per_sequence.tsv")
os.system("sort -k2 " + outdir + "taxonomy_assignment_per_sequence_eight_levels.tsv > " + outdir + "tmp.txt && mv " + outdir + "tmp.txt " + outdir + "taxonomy_assignment_per_sequence_eight_levels.tsv")



#Done
log_and_print('Taxonomy assignment Complete. Have a nice day!')
