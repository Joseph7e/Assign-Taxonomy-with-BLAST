# Assign-Taxonomy-with-BLAST
This script can be used for any standard multi-fasta.
Will work on any database, including custom ones or ncbi's nt database.
A path to taxonomy classifications is required but several precomputed are available here.


# Dependencies
###PYTHON3
modules: Biopython
### For commputing otus
pick_otus.py

http://qiime.org/scripts/pick_otus.html
### Databases
ncbi's nt database:

SILVA 18S databases: https://www.arb-silva.de/download/archive/qiime/
### Taxonomy classifications
ncbi taxa dump: expanded and customized for this program: AVAILABLE ON THIS GIT FOR DOWNLOAD

SILVAs taxonomy databases are available at SILVA link above

# INSTALLATION
Download the script from above (right click above link, copy link address)
wget 


# How It Works

### OTU PICKING
The first step is OTU picking. This is done to limit the number of time consuming blast hits that need to be done.
The default is to use qiimes otu_picking.py script which provides the correctly formatted otu_seqs.txt file.
A seed sequence is chosen from the OTU text file and used to represent that OTU (currently selecting the largest sequence belonging to that OTU cluster). Those seeds are then the only sequences that are blasted.
This can change the number of blasts down drastically. The original sequences can then be taxified based on the OTU seeds blast result.
I WILL ADD AN OPTION TO AVOID OTU PICKING IN THE FUTURE IF REQUESTED.

### Taxonomy Assignment with BLAST

#### Consensus assignment
Rather than relying on a single blast hit the program takes the top X (user defined with --hits_to_consider) blast hits for each sequence.
It then computes the classification for each of these blast hits based on the best consensus taxonomy. If all ten blast hits agree on the same taxonomy than it will give the value to the species level, however if they only agree to the family level then it will stop there.

This aspect has many options, including setting the maximum number of blast hits to consider and the percent sway from the best blast hit. Blast hits are not all treated the same. If your blast provided 10 blast hits but only 5 of them are within 'percent_sway' the others will not be considered. For example, if one blast hit has a percent identity of 99% while the others are only 90%, only the top hit will be considered (unless of course you set the --percent_sway option to 9.0 or above). 

In addition the program will mask blast hits that have a taxonomic assignment of uncultured, or unclassified. These hits will only be considered as a last resort.

If you want to only consider the best blast hit, just set --hits_to_consider to '1'.

#### Best taxonomy based on percent identity
The program provides another level of taxonomic certainty based on the blast percentage. For example, if the best blast hit is 90% you can be fairly confident that you cannot provide the species of the organisms but maybe you can provide the phylum. Right now there are three levels that you can set with the program, --cutoff_species, --cutoff_family, and --cutoff_phylum. The phylum level cutoff is also used as an ultimate filter for the blast hits.


#### DOING YOUR OWN BLAST?? We will import that and save you the step.
use this  -outfmt '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids'
MAKE SURE THAT YOUR BLAST IS RUN WITH THE REP_OTU_SEQS.FASTA, I WILL ADD OPTION TO USE PRE-COMPUTED SEQS.FNA BLAST IN THE FUTURE.


