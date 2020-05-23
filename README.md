# Project Title
DNCON

## What it does?
This project is reconstruction of **DNCON:Predicting protein residue–residue contacts using deep networks and boosting**. The base paper link is here: https://academic.oup.com/bioinformatics/article/28/23/3066/195693#92186201.

## Extra notes
Every files test cases should be bottom of that file.

## Directories, modules and packages
```
configs/
    __init_.py
    general_config.py
data/
    contact_maps/
    fastas/
    pdbs/
datasets/
    __init_.py
metrics/
    __init_.py
modules/
    __init_.py
outputs/
    images/
    logs/
    models/
utils
    __init_.py
vizualization
    __init_.py
    data_viz.py
    output_viz.py
run.ipnyb
run.py
tester.ipnyb
tester.py
```
## Requirements
Python 3

## How to run?
When you first use this, try running the following command to see if all packages and modules are successfully imported.
```
>> python tester.py
Successfully imported all packages
```
You could use the following command:
```
>> python run.py
```
Or you can open jupyter notebook to run the notebooks.

## Todo
1. **PSSM using PSI-BLAST**<br />
    First BLAST command line applications (i.e. psiblast, blastp etc) are used using the following manner. Then ```Bio.Blast.Applications``` is used for running ```psiblast``` using Biopython.
    <br />
    Quick start on BLAST: https://www.ncbi.nlm.nih.gov/books/NBK279680/
    <br />
    BLAST TOOLs: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ <br />
    BLAST DBs: https://ftp.ncbi.nlm.nih.gov/blast/db/ <br />
    Currently the work is done using "swissprot" database. Later it must be done using  a non-redundant version of the nr sequence database.
    <br />

    To download "swissprot":
    ```
    ncbi/ncbi-blast-2.10.0+/bin/update_blastdb.pl --decompress swissprot
    ```
    then keep all files in "/my_swissprot/" directory

    To download all "nr" databases: 
    ```
    update_blastdb.pl --decompress nr [*]
    ```

    To download "swissprot" fasta file go to: https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/ <br />
    To preformat downloaded "swissprot.tar.gz" go to: https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/ <br />

    The following command can be used to create blast-db from fasta sequence
    ```
    ncbi/ncbi-blast-2.10.0+/bin/makeblastdb -dbtype prot -in swissprot.fasta -input_type fasta -out swissprot_mine/swissprot

    Here, swissprot is the name of the database (not the path), couple of files names are prefixed with this.
    ```

2. **Scratch Protein predictor** for secondary structure (SSpro) and solvent accessibility (ACCpro) computation: http://scratch.proteomics.ics.uci.edu/
3.


## Notes

