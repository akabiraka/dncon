# data directories
DATA_DIR = "data/"
ALL_PDB_IDS = DATA_DIR + "all_pdb_ids_tiny.txt"
PDB_DIR = DATA_DIR + "pdbs/"
BLAST_DIR = DATA_DIR + "blast/"
PSIBLAST_DIR = DATA_DIR + "psiblast/"
CONTACT_MAP_DIR = DATA_DIR + "contact_maps/"
FASTA_DIR = DATA_DIR + "fastas/"
BAD_PDB_IDS = DATA_DIR + "bad_pdb_ids.txt"
GOOD_PDB_IDS = DATA_DIR + "good_pdb_ids.txt"

# BLAST directories
NCBI = "ncbi/"
BLAST_DIR = NCBI + "ncbi-blast-2.10.0+/"
PSIBLAST_EXE = BLAST_DIR + "bin/psiblast"
BLAST_DB = NCBI + "my_swissprot/swissprot"

# file format
CIF = 'mmCif'
CIF_ATOM = "cif-atom"
FASTA = 'fasta'

# file extensions
DOT_CIF = ".cif"
DOT_PT = ".pt"
DOT_FASTA = ".fasta"
DOT_XML = ".xml"

# output directories
OUTPUT_DIR = "outputs/"
OUTPUT_IMAGES_DIR = OUTPUT_DIR + "images/"
OUTPUT_MODELS_DIR = OUTPUT_DIR + "models/"
OUTPUT_LOGS_DIR = OUTPUT_DIR + "logs/"

# amino-acids
AMINO_ACID_3TO1 = {'ALA': 'A',
                   'CYS': 'C',
                   'ASP': 'D',
                   'GLU': 'E',
                   'PHE': 'F',
                   'GLY': 'G',
                   'HIS': 'H',
                   'ILE': 'I',
                   'LYS': 'K',
                   'LEU': 'L',
                   'MET': 'M',
                   'ASN': 'N',
                   'PRO': 'P',
                   'GLN': 'Q',
                   'ARG': 'R',
                   'SER': 'S',
                   'THR': 'T',
                   'VAL': 'V',
                   'TRP': 'W',
                   'TYR': 'Y'}

