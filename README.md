# Assign-Taxonomy-with-BLAST
This program was primarily designed for the taxonomic assignment of amplicon sequence variants (ASVs), but it works with any sequence data in FASTA format.
Many output files are generated, several of which can be directly imported into qiime or other common workflows.

# Dependencies
###Python3
tested with Python 3.4.3
modules: Biopython, argsparse

## Compute ASVs/OTUs (optional)
If your starting data is fastq data from an amplicon based experiment you will need to compute ASVs/OTUs prior to running this program. Here are some great options.
dada2 - https://benjjneb.github.io/dada2/tutorial.html
qiime2 -  
mothur - 

## Prepare sequence database and taxonomic lookup
ncbi's nt database - ftp://ftp.ncbi.nlm.nih.gov/blast/db/
```bash
#####To download this database copy and paste the line below (it will take a bit)
mkdir ncbi_nt_database && cd ncbi_nt_database
wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.gz'
tar -xvf nt*
```

SILVA ssuRNA sequence and taxonomy database: https://www.arb-silva.de/download/archive/qiime/

##### Constructing an updated taxonomy database for use with NCBI nt.

Locate and download the latest NCBI taxonomy database. 
You will need the names and nodes file to construct the expanded taxonomy lookup.
Link updated December, 2020

```
mkdir ncbi_taxonomy/ && cd ncbi_taxonomy/
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xvzf *.tar.gz
```

Expand the taxonomic lineages into a simplified tsv lookup.
```
python3 genbank_nodes_and_names_to_taxonomy.py  names.dmp nodes.dmp
```


# How It Works

#### Consensus assignment
Rather than relying on a single blast hit the program takes the top X (user defined with --hits_to_consider) blast hits for each sequence.
It then computes the classification for each of these blast hits and determines taxonomy for the query sequence based on the best consensus taxonomy. i.e. If all ten blast hits agree on the same taxonomy than it will give the value to the species level, however if they only agree to the family level then it will stop there.

This aspect has many options, including setting the maximum number of blast hits to consider and the percent sway from the best blast hit. Blast hits are not all treated the same. If your blast provided 10 blast hits but only 5 of them are within '--percent_sway' the others will not be considered. For example, if one blast hit has a percent identity of 99% while the others are only 95%, only the top hit will be considered (unless of course you set the --percent_sway option to 4.0 or above). The default --percent_sway value is 0.5, this will gather very similiar hits and ignore those that are less similiar.

If you want to only consider the best blast hit, just set --hits_to_consider to '1'.


#### Best taxonomy based on percent identity
The program provides another level of taxonomic certainty based on the blast percentage. For example, if the best blast hit is 90% you can be fairly confident that you cannot provide the species of the organisms but maybe you can provide the phylum. Right now there are three levels that you can set with the program, --cutoff_species, --cutoff_family, and --cutoff_phylum. The phylum level cutoff is also used as an ultimate filter for the blast hits percent identity.

If you want to only keep sequences that are identified to the species level, just set all cutoffs to '97' or '99'.

If you want to leave it to the consensus taxonomy to decide 'best estimated taxonomy', set all three cutoffs low, maybe '80'.

Note: The consensus taxonomy usually does a pretty good job of weeding out incorrect taxonomy, setting all these cutoff values to 80 actually provides pretty great taxonomy in a lot of cases. Especially if you set a high --hits_to_consider and --percent_sway.

#### Masking uncultured taxonomy
By default the program will mask blast hits that have a taxonomic assignment of uncultured, or unclassified. These hits will only be considered as a last resort.


#### Other blast options
You have the option to set the minimum length coverage for the blast hit (defined by 'length of query'/'length of hit'). The default is 0.8. I like to be stringent here to avoid tiny insignificant blast hits.


#### You already did YOUR OWN BLAST?? We will import that and save you the step. Just make sure the format is right.
use this format --> -outfmt '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids'
If running your own blast command make sure you run it with the out_seqs.fasta file. I have not added an option to input seqs.fna blast directly into qiime yet… but will if requested.


#USAGE
python3 taxonomy_assignment_BLAST.py [options] sequence_file blast_database taxonomy_file  
### To see all options
```python3 taxonomy_assignment_BLAST.py -h  

usage: taxonomy_assignment_BLAST_V2.py [-h] [-v]  
                                       [--cutoff_species CUTOFF_SPECIES]  
                                       [--cutoff_family CUTOFF_FAMILY]  
                                       [--cutoff_phylum CUTOFF_PHYLUM]  
                                       [--length_percentage LENGTH_PERCENTAGE]  
                                       [--length_cutoff LENGTH_CUTOFF]  
                                       [--hits_to_consider HITS_TO_CONSIDER]  
                                       [--percent_sway PERCENT_SWAY]  
                                       [--blast_evalue BLAST_EVALUE]  
                                       [--blast_threads BLAST_THREADS]  
                                       [--blast_flavor BLAST_FLAVOR]  
                                       [--blast_database BLAST_DATABASE]  
                                       [--blast_file BLAST_FILE] [--ncbi_nt]  
                                       [--output_dir OUTPUT_DIR]  
                                       [--config_file CONFIG_FILE]  
                                       sequence_file tax_file  

positional arguments:  
  sequence_file         seqs.fna file from qiime or any multifasta, just make  
                        sure header has unique id with a space  
  tax_file              path to silva or customized blast database taxonomy  
                        file  
  
optional arguments:  
  -h, --help            show this help message and exit  
  -v, --verbose         increase output verbosity (default: False)  
  --cutoff_species CUTOFF_SPECIES  
                        cutoff for finest taxonomic level (default: 97)  
  --cutoff_family CUTOFF_FAMILY  
                        cutoff for family taxonomic level (default: 90)  
  --cutoff_phylum CUTOFF_PHYLUM  
                        cutoff for phylum taxonomic level, also acts as  
                        ultimate cutoff value for blast (default: 80)  
  --length_percentage LENGTH_PERCENTAGE  
                        cutoff for query_hit/length_of_query i.e. query  
                        coverage (default: 0.0)  
  --length_cutoff LENGTH_CUTOFF  
                        primary cutoff for length of hit (default: 0)  
  --hits_to_consider HITS_TO_CONSIDER  
                        number of hits to consider when gathering consensus  
                        taxonomy (default: 3)  
  --percent_sway PERCENT_SWAY  
                        when comparing greater than 1 blast hit, value is the  
                        percent from best hit considered. i.e. if best blast  
                        hit is 99.5 ID, a value of 0.5 will consider  
                        everything 99 and greater when creating the consensus  
                        (default: 0.5)  
  --blast_evalue BLAST_EVALUE  
                        setting for e-value cutoff for blast, must be in form  
                        1e-X (default: 1e-10)  
  --blast_threads BLAST_THREADS  
                        set the number of threads for blast (default: 24)  
  --blast_flavor BLAST_FLAVOR  
                        select blastp, tblastn etc, must have correct database  
                        and query sequence formats (default: blastn)  
  --blast_database BLAST_DATABASE  
                        path tp blast database, be sure to run makeblastdb on  
                        the database first, type 'IGNORE' if precomputed blast  
                        is given (default: None)  
  --blast_file BLAST_FILE  
                        precomputed blast results, MUST BE MY CUSTOMIZED  
                        FORMAT, '6 qseqid qlen sseqid pident length qstart  
                        qend sstart send evalue bitscore' (default: None)  
  --ncbi_nt             REQUIRED flag for use of ncbi nt database (default:False)  
  --output_dir OUTPUT_DIR  
                        output directory name (default: Assigned_Taxonomy)  
  --config_file CONFIG_FILE  
                        NOT YET IMPLEMENTED, an optional parameters file, any  
                        arguments in this file will override other options  
                        (default: None)  ```
