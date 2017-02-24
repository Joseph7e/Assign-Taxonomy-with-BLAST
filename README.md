# Assign-Taxonomy-with-BLAST
This script is used for otu clustering, blasting, and ultimately, assigning taxonomy. All steps are optional and customizable.

It can be used for any multi-fasta and should work on any database, including custom ones or ncbi's nt database.

A path to taxonomy classifications for the database is required but several have already been precomputed and are available here. If you would like me to construct a taxonomy classification for your database or instructions on how to do so. Email me.

Many output files are generated. Several of which can be directly imported into qiime or qiime2 for visualization of taxonomy or other analyses.




# Dependencies
###Python3
tested with Python 3.4.3

modules: Biopython
### Programs for computing otus (currently required)
pick_otus.py
available here: --> http://qiime.org/scripts/pick_otus.html
#### other options coming soon

### Databases (not required)
ncbi's nt database: --> ftp://ftp.ncbi.nlm.nih.gov/blast/db/
####To download this database
mkdir ncbi_nt_database && cd ncbi_nt_database && wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.gz' && tar -xvf nt*

SILVA 18S databases: https://www.arb-silva.de/download/archive/qiime/
### Taxonomy classifications (you will need one for your database)
ncbi taxa dump: expanded and customized for this program: AVAILABLE ON THIS GIT FOR DOWNLOAD

SILVAs taxonomy databases: --> https://www.arb-silva.de/download/archive/qiime/

# Installation
Download the script from above. Thats it.

If you want the ncbi precomputed taxonomic file download from --> http://cobb.unh.edu/ncbi_taxonomy_expanded.tsv.gz

Then gunzip the file.

# How It Works

### OTU PICKING
The first step is OTU picking. This is done to limit the number of time consuming blasts that need to be done, and for prepping the sequences for qiime.

The default is to use uclust via the otu_picking.py script, which provides the correctly formatted otu_seqs.txt file.

Next the program will choose a seed sequence from the resulting OTU text file info, this seed will be used to represent that OTU (currently selecting the largest sequence belonging to that OTU cluster) and is the only sequence that is blasted for that cluster.
This often cuts the number of blasts down drastically. The original sequences can then be taxified based on the OTU seeds blast result.
I can add an option to ignore otu picking if requested.

### Taxonomy Assignment with BLAST

#### Consensus assignment
Rather than relying on a single blast hit the program takes the top X (user defined with --hits_to_consider) blast hits for each sequence.
It then computes the classification for each of these blast hits and determines taxonomy for the OTU based on the best consensus taxonomy. i.e. If all ten blast hits agree on the same taxonomy than it will give the value to the species level, however if they only agree to the family level then it will stop there.

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
python3 taxonomy_assignment_BLAST.py -h

### Info about required input formats

#### sequence_file
#### blast_database
#### taxonomy_file

### Info about optional input files

#### blast_file
#### otu_file

### Info about output
