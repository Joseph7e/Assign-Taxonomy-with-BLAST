#!/usr/bin/python3
import sys

names_file = sys.argv[1] #'/mnt/lustre/hcgs/joseph7e/databases/ncbi_taxonomy_db/names.dmp'#sys.argv[1]
nodes_file = sys.argv[2] #'/mnt/lustre/hcgs/joseph7e/databases/ncbi_taxonomy_db/nodes.dmp' #sys.argv[2]

names_dictionary = {} # tax_id: name_text
nodes_dictionary = {} # tax_id: rank, parent_tax_id

#master_dict[id] = [rank,name_text,parent_tax_id]

print ('Parsing names.dmp file for tax_id and scientific_name...')
for l in open(names_file):
    tax_id, name_text, unique_name, name_class, tmp = [x.lstrip().rstrip() for x in l.rstrip().split('|')]
    if tax_id not in names_dictionary.keys() and name_class == 'scientific name' : # only get one item added for each tax_id using valid id (scientific name)
        names_dictionary[tax_id] = name_text
        #print (tax_id, name_text, unique_name, name_class)
print ('A total of {} unique entries in the names database'.format(len(names_dictionary.keys())))
print ("    (this number should equal -->   awk -F'\t|' '{print $1}' names.dmp | uniq | wc        )")


print ('Parsing nodes.dmp for tax_id, parent_tax_id and taxonomy rank...')
print ('    (this file also has genetic code, mito_gen_cod and parent codes as well)')
for n in open(nodes_file):
    node_elements = [x.lstrip().rstrip() for x in n.rstrip().split('|')]
    #tax_id, parent_tax_id, rank, embl_code, division_id, tmp, gen_code, tmp, mito_code, tmp, tmp, tmp, comments, tmp = [x.lstrip().rstrip() for x in n.rstrip().split('|')]
    tax_id, parent_tax_id, rank = node_elements[0:3]
    if tax_id not in nodes_dictionary.keys():
        nodes_dictionary[tax_id] = [rank, parent_tax_id]
print ('A total of {} unique entries in the names database'.format(len(nodes_dictionary.keys())))


def Expand_taxonomy(tax_id, names_dictionary, nodes_dictionary):
    tax_levels = ['superkingdom', 'kingdom', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'superorder', 'order', 'superfamily', 'family', 'subfamily', 'genus', 'species']
    growing_taxonomy = [names_dictionary[tax_id]]
    growing_taxonomy_levels = [nodes_dictionary[tax_id][0]]
    master_key = tax_id

    while True:
        new_tax_id = nodes_dictionary[tax_id][1] # parent_tax_id
        new_name = names_dictionary[new_tax_id]
        new_rank = nodes_dictionary[new_tax_id][0]
        if new_name == 'root':
            break

        growing_taxonomy.append(new_name)
        growing_taxonomy_levels.append(new_rank)
        tax_id = new_tax_id
    reversed_taxonomy = list(reversed(growing_taxonomy))
    reversed_taxonomy_levels = list(reversed(growing_taxonomy_levels))

    #final_tax = ['unknown_superkingdom', 'unknown_kingdom', 'unknown_phylum', 'unknown_subphylum', 'unknown_superclass', 'unknown_class', 'unknown_subclass', 'unknown_superorder', 'unknown_order', 'unknown_superfamily', 'unknown_family', 'unknown_subfamily', 'unknown_family', 'unknown_genus', 'unknown_species']
    final_tax = ['unknown_' + x for x in tax_levels]

    for t in range(len(tax_levels)):
        for i in range(len(reversed_taxonomy_levels)):
            if reversed_taxonomy_levels[i] == tax_levels[t]:
                final_tax[t] = reversed_taxonomy[i]

    # format_tax = []
    # for i in range(len(final_tax)):
    #     format_tax.append('D_'+str(i)+'__'+final_tax[i])
    #combined_taxonomy = ';'.join(format_tax)

    combined_taxonomy = ';'.join(final_tax)
    #print (master_key+'\t'+combined_taxonomy)
    return master_key, combined_taxonomy




output_file = open('expanded_ncbi_taxonomy.tsv', 'w')

for tax_id in list(names_dictionary.keys()):#.sort(key=int):
    t, expanded_taxonomy = Expand_taxonomy(tax_id, names_dictionary, nodes_dictionary)
    output_file.writelines(t + '\t' + expanded_taxonomy+'\n')


