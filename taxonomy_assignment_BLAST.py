#!/usr/bin/python3
#Author: Joseph Sevigny
#Affiliation: Hubbard Center for Genome Studies, University of New Hampshire
#Date: 02/23/2017

import sys, re, os, subprocess, argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#OPTIONAL ARGUMENTS
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
#parser.add_argument("-h","--help","-help", help="prints help and exits")
parser.add_argument("--cutoff_species", help="cutoff for finest taxonomic level", type=int, default=97)
parser.add_argument("--cutoff_family", help="cutoff for family taxonomic level", type=int, default=90)
parser.add_argument("--cutoff_phylum", help="cutoff for phylum taxonomic level, also acts as ultimate cutoff value for blast", type=int, default=80)
parser.add_argument("--length_percentage", help="cutoff for query_hit/length_of_query i.e. query coverage", type=float, default=0.8)
parser.add_argument("--length_cutoff", help="primary cutoff for length of hit", type=int, default=0)
parser.add_argument("--hits_to_consider", help="number of hits to consider when gathering consensus taxonomy", type=float, default=5)
parser.add_argument("--percent_sway", help="when comparing greater than 1 blast hit, value is the percent from best hit considered. i.e. if best blast hit is 99.5 ID, a value of 0.5 will consider everything 99 and greater when creating the consensus", type=float, default=0.5)
parser.add_argument("--otu_file", help="precomputed otu file, will run otu picking if not given")
parser.add_argument("--blast_file", help="precomputed blast results, MUST BE MY CUSTOMIZED FORMAT, '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore'")
parser.add_argument("--ncbi_nt", help="flag for use of ncbi nt database", action="store_true")
parser.add_argument("--blast_evalue", help="setting for e-value cutoff for blast, must be in form 1e-X", type=str, default="1e-10")
parser.add_argument("--make_biom", help="make a otu biom table using otus and taoxnomy assignments", action="store_true")

#REQUIRED ARGUMENTS
parser.add_argument("sequence_file", help="seqs.fna file from qiime or any multifasta, just make sure header has unique id with a space")
parser.add_argument("blast_database", help="path tp blast database, be sure to run makeblastdb on the database first, type 'IGNORE' if precomputed blast is given")
parser.add_argument("tax_file", help="path to silva or customized blast database taxonomy file")


#parser.add_argument
args = parser.parse_args()
if args.verbose:
    print ('####VERBOSE OPTION GIVEN\n')
    print ('MAIN ARGUMENTS')
    print ('sequence_file: ', args.sequence_file)
    print ('blast_database: ', args.blast_database)
    print ('taxonomy_file: ', args.tax_file)
    print ('\nOPTIONS')
    print ('blast_hits_to_consider: ', args.hits_to_consider)
    print ('sway_from_best_hit: ', args.percent_sway)
    print ('phylum_level_cutoff: ', args.cutoff_phylum)
    print ('family_level_cutoff: ', args.cutoff_family)
    print ('species_level_cutoff', args.cutoff_species)
    print ('length_of_hit_cutoff: ', args.length_cutoff)
    print ('proportion_of_query_coverage: ', args.length_percentage)
    print ('blast e-value cutoff: ', args.blast_evalue)
    print ('Construct OTU BIOM table: ', args.make_biom)
    print ('OTU_file given?: ', args.otu_file)
    print ('BLAST_file given?: ', args.blast_file)
    print ('NCBI flag?: ', args.ncbi_nt)


###INFO FOR TAXONOMY DATABASES
taxonomy_categories = ['superkingdom', 'kingdom', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'superorder', 'order', 'superfamily', 'family', 'subfamily', 'family', 'genus', 'species']
#qiime_all_level_12_categories = [superkingdom, subkingdom, sub_subkingdom, kingdom, tmp1, tmp2, phylum, subphylum, class, order, family, genus, species]

### OTHER OPTIONS BASED ON DATABASE OF CHOICE
uncultured_cutoff_level = len(taxonomy_categories)
species_level = 14
family_level = 10
phylum_level = 2
#taxonomy_to_grab = [0,3,6,10,11,12] # file for qiime output
taxonomy_to_grab = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]

#OUTPUT_FILES
outdir = 'Assigned_Taxonomy/'; os.mkdir(outdir)
otu_sequences_output = outdir + 'otu_seqs.fasta'
taxonomy_assignment_outfile = outdir + 'taxonomy_assignment_per_sequence.txt'
log_file_handle = open(outdir + 'log_file.txt','w')



#Potential Additions:
## --> Add option to blast more than one otu if they might be different, make the otu selection strigent




#####Dictionaries

sequence_dict = SeqIO.to_dict(SeqIO.parse(args.sequence_file, "fasta"))
otu_dictionary = {} # otu:[sequences]
taxonomy_dictionary = {} # tax_code:[taxonomy]
otu_taxonomy_dict = {} # otu:[best_taxonomy]
percent_id_dict = {}


######OTU PICKING

def run_otu_picker(sequence_file, otu_output_dir):
    '''perform otu_picking with uclust'''
    otu_pid = '0.99'
    otu_command = ["pick_otus.py", "-i", sequence_file, "-o", otu_output_dir, "-s", otu_pid, "--enable_rev_strand_match"] #command to complete otus
    # if args.verbose:
    #     print ('Command =', ' '.join(otu_command))
    sp = subprocess.Popen(otu_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)#blast = sp.communicate()
    sp.communicate()



if args.otu_file == None:
    if args.verbose:
        print ('\nConstructing OTUS from sequence file')
    otu_output_dir = "otu_picking_output/"
    run_otu_picker(args.sequence_file, otu_output_dir)
    otu_list_dir = os.listdir(otu_output_dir)
    for file in otu_list_dir:
        if file.endswith('.txt'):
            otu_file = otu_output_dir + file
else:
    otu_file = args.otu_file


###parse otu file

if args.verbose:
    print ('Constructing OTU seed fasta file for blast')
for line in open(otu_file,'r'):
    elements = line.rstrip().split('\t')
    group = elements[0]
    for sample_id in elements[1:]:
        if group in otu_dictionary.keys():
            otu_dictionary[group].append(sample_id)
        else:
            otu_dictionary[group] = [sample_id]

###parse sequence_file and pick otu seed (select reference fasta for blast)
otu_sequences_output_handle = open(outdir + 'otu_seqs.fasta', 'w')
for key, value in otu_dictionary.items(): #could make this a little better than first one????
    otu_sequences_output_handle.writelines('>'+key+'\n' + sequence_dict[value[0]].seq+'\n')



#parse taxonomy file

if args.verbose:
    print ('Parsing taxonomy file and storing info in dictionary (might take a moment if using nt database)')
for line in open(args.tax_file,'r'):
    elements = line.rstrip().split('\t')
    tax_id = elements[0]
    tax_info =re.sub('D_[0-9]*__','',elements[1]).replace(' ','_').split(';')
    taxonomy_dictionary[tax_id]=tax_info
    #taxonomy_dictionary[tax_id]=[tax_info[3],tax_info[6:]] ???

###do_blast
def do_blast(query, database, blast_output):
    '''performs blasts and returns the results'''
    custom_blast_format = '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids'
    #db_command = ["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", "temp_db"] #command to construct database
    #sd = subprocess.Popen(db_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE); sd.communicate()
    blast_command = ["blastn", "-query", query, "-db", database , "-num_threads", "48", "-evalue", args.blast_evalue, "-out", blast_output, "-outfmt", custom_blast_format] #command to complete blastn
    sp = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)#blast = sp.communicate()
    sp.communicate()

if args.blast_file == None:
    if args.verbose:
        print ('Running blast on OTU seeds (this step takes the longest, maybe run overnight')
    blast_file = 'otu_seqs_blast.tsv'
    do_blast(otu_sequences_output,args.blast_database,blast_file)
else:
    blast_file = args.blast_file


####PARSE BLAST and ASSIGN TAXONOMY#######
#make so assignment avoids one mis when the rest agree

if args.verbose:
    print ('Parsing Blast and Assigning Taxonomy')

current_best_hits = {} #count_code:[percent_id, tax_code]
current_query = ''
for line in open(blast_file,'r'):
    elements = line.rstrip().split('\t'); sequence_id = elements[0]; tax_code = elements[2]
    percent_id = elements[3]; length_hit = elements[4]; length_query = elements[1]; ncbi_id = elements[-1]

    if args.ncbi_nt:
        if ';' in ncbi_id:
            ncbi_id = ncbi_id.split(';')[0]
        tax_code = ncbi_id

    if current_query != sequence_id: #INITIATION of new query
        if bool(current_best_hits) and current_query: #ENSURE A SET OF BLASTS TO PARSE

            if args.verbose:
                print ('\nASSIGNING TAXONOMY FOR', current_query, 'total hits passing initial filters = ', len(current_best_hits))

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
                    if args.verbose:
                        print (current_query, value[1], value[0], '--> CAPTURED after percent sway filter')
            ##Assign Taxonomic Certainty
            level = phylum_level #'phylum'
            if max(top_hits) >= args.cutoff_family:
                level = family_level #'family'
            if max(top_hits) >= args.cutoff_species:
                level = species_level #'species'



            blast_percent = max(top_hits)

            ##Assign Taxonomy
            temp_taxonomy_dict = {} #unique_code:[taxonomy]
            count = 0
            for key,value in best_blast_dictionary.items():
                count += 1
                pkey = str(value[0])+ '_'+str(count)
                #Mask uncultured taxonomy
                if 'uncultured' in ' '.join(taxonomy_dictionary[value[1]][:uncultured_cutoff_level+1]) or 'unidentified' in ' '.join(taxonomy_dictionary[value[1]][:uncultured_cutoff_level+1]):
                    pkey = 'X'+pkey
                    #print ('UNCULTURED', current_query)
                temp_taxonomy_dict[pkey]=taxonomy_dictionary[value[1]]#[:level]


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
                    best_level_taxonomy = value[:level]

            best_level_taxonomy = list(filter(bool, best_level_taxonomy))#remove empty taxonomic levels
            needed_unknowns = (len(taxonomy_categories)-len(best_level_taxonomy))
            for i in range(needed_unknowns):
                best_level_taxonomy.append('undetermined')

            ### FINALLY FILL OTU INFORMATION!!!!
            otu_taxonomy_dict[current_query] = best_level_taxonomy
            percent_id_dict[current_query] = blast_percent
            if args.verbose:
                print ('Taxonomy Assignment for', current_query, '=', best_level_taxonomy[0], best_level_taxonomy[phylum_level], best_level_taxonomy[family_level], best_level_taxonomy[-1])

        #reset query_id
        current_query = sequence_id
        current_best_hits = {}

    ## Filter Hits and Loop Through current Query and add up to x number of blast hits
    if int(length_hit)/int(length_query)>args.length_percentage and float(percent_id)>args.cutoff_phylum and int(length_hit)>args.length_cutoff:
        if len(current_best_hits.keys()) < args.hits_to_consider:
            label = len(current_best_hits.keys()) # provide unique key based on number of matched
            current_best_hits[label] = [percent_id, tax_code]



### QIIME FORMATTED OUTPUT

if args.verbose:
    print ('\nConstructing output files')

new_output_file = outdir + 'taxonomy_assignment_qiime_format.txt'

for otu, sequences in otu_dictionary.items():
    output_line = otu+'\t'
    if otu in otu_taxonomy_dict.keys():
        blast_percent = round((percent_id_dict[otu]/100),2)
        count = 0
        for taxonomy_level in otu_taxonomy_dict[otu]:
            if count in taxonomy_to_grab:
                output_line += 'D_'+str(count)+'__'+taxonomy_level+';'
                count +=1
        output_line = output_line[:-1] + '\t'+str(blast_percent)+'\t1'
    else:
        output_line += 'Unassigned\t1.00\t1'

    with open(new_output_file,'a') as t:
        t.writelines(output_line+'\n')


#Regular Assignment

for otu, sequences in otu_dictionary.items():
    for seq in sequences:
        if '__' in seq:
            group,seq_id = seq.split('__')
        else:
            group,seq_id = seq.split('_') #if two of these _ are in it
        output_line = group+':'+seq_id+':'+otu
        if otu in otu_taxonomy_dict.keys():
            for taxonomy_level in otu_taxonomy_dict[otu]:
                output_line += (':' + taxonomy_level)
            with open(taxonomy_assignment_outfile,'a') as t:
                t.writelines(output_line+'\n')

        else:
            output_line += (':undetermined'*len(taxonomy_to_grab))
            with open(taxonomy_assignment_outfile,'a') as t:
                t.writelines(output_line+'\n')

## sort the taxonomy file for viewing
os.system("sort -k2 Assigned_Taxonomy/taxonomy_assignment_qiime_format.txt > Assigned_Taxonomy/taxonomy_assignment_qiime_format2.txt && mv Assigned_Taxonomy/taxonomy_assignment_qiime_format2.txt Assigned_Taxonomy/taxonomy_assignment_qiime_format.txt")
os.system("sort -t: -k4 Assigned_Taxonomy/taxonomy_assignment_per_sequence.txt > Assigned_Taxonomy/taxonomy_assignment2.txt && mv Assigned_Taxonomy/taxonomy_assignment2.txt Assigned_Taxonomy/taxonomy_assignment_per_sequence.txt")


#Targeted Assignment

# for otu, sequences in otu_dictionary.items():
#     #print (sequences)
#     for seq in sequences:
#         if '__' in seq:
#             group,seq_id = seq.split('__')
#         else:
#             group,seq_id = seq.split('_') #if two of these _ are in it
#         output_line = group+':'+seq_id+':'+otu
#         if otu in otu_taxonomy_dict.keys():
#             #print (otu_taxonomy_dict[otu][1])
#             if otu_taxonomy_dict[otu][3] == 'Metazoa_(Animalia)':
#                 for taxonomy_level in otu_taxonomy_dict[otu]:
#                     output_line += (':' + taxonomy_level)
#
#             elif otu_taxonomy_dict[otu][3] == 'Fungi':
#                 output_line += (':Fungi'*12)
#             else:
#                 #BREAK INTO Taxonomic Groups
#                 non_metazoa = ':non-metazoa'
#                 output_line += ((':'+otu_taxonomy_dict[otu][1])+non_metazoa*11)
#
#         else:
#             output_line += (':undetermined'*12)
#
#         with open(taxonomy_assignment_outfile,'a') as t:
#             t.writelines(output_line+'\n')





def make_otu_table(otu_file,taxonomy_file):
    '''makes otu_table for qiime'''
    output = 'otu_table.biom'
    biom_command = ["make_otu_table.py", "-i", otu_file, "-t", taxonomy_file, "-o", output ] #command to complete blastn
    sp = subprocess.Popen(biom_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)#blast = sp.communicate()
    sp.communicate()
    os.system('mv otu_table.biom Assigned_Taxonomy/')

if args.make_biom:
    if args.verbose:
        print ('Constructing OTU table for qiime')
    make_otu_table(otu_file, new_output_file)

if args.verbose:
    print ('Process Complete. Have a nice day!')
