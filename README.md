# Assign-Taxonomy-with-BLAST
This script can be used for any standard multi-fasta.
Will work on any database, including custom ones or ncbi's nt database.
A path to taxonomy classifications is required but several precomputed are available here.


# Dependencies
not much
### PYTHON
python3
modules: Biopython

### For commputing otus
pick_otus.py
http://qiime.org/scripts/pick_otus.html

### Databases
ncbi's nt database: \n
SILVA 18S databases: https://www.arb-silva.de/download/archive/qiime/

### Taxonomy classifications
ncbi taxa dump: expanded and customized for this program: AVAILABLE ON THIS GIT PAGE\n
SILVAs taxonomy databases are available at SILVA link above


## INSTALLATION
Just download the script, thats it.


# How It Works
### DOING YOUR OWN BLAST??
use this as -outfmt '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids'


