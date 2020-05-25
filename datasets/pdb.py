import sys
sys.path.append('../dncon')
import numpy as np
import traceback
import subprocess

from Bio.PDB import *
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline

import configs.general_config as CONFIGS
import utils.data_utils as DataUtils
import vizualizations.data_viz as DataViz


class PDB(object):
    def __init__(self):
        print("alhumdulillah")
        self.pdb_ids_file = CONFIGS.ALL_PDB_IDS

        self.pdbl = PDBList()
        self.parser = MMCIFParser(QUIET=True)
        self.aa_3to1 = CONFIGS.AMINO_ACID_3TO1

        self.threshhold = 8

        self.blast_db = CONFIGS.BLAST_DB
        self.psiblast_exe = CONFIGS.PSIBLAST_EXE
        
    def read_pdb_ids(self, line):
        """
            Assuming a line has 6 columns in a row.
            returns:
                pdb_code, chain_id: ('5SY8', 'O') in upper case
        """
        line = line.split()
        pdb_code = line[0][:-1].lower()
        chain_id = line[0][-1]
        return pdb_code, chain_id

    def download(self, pdb_code):
        # downloading in .cif format
        self.pdbl.retrieve_pdb_file(pdb_code, pdir=CONFIGS.PDB_DIR, file_format=CONFIGS.CIF)

    def filter_aa_residues(self, chain):
        """
        a chain can be heteroatoms(water, ions, etc; anything that isn't an amino acid or nucleic acid)
        so this function get rid of atoms excepts amino-acids
        """
        aa_residues = []
        non_aa_residues = []
        non_aa = []
        seq = ""
        for i in chain:
            if i.get_resname() in standard_aa_names:
                aa_residues.append(i)
                seq += self.aa_3to1[i.get_resname()]
            else:
                non_aa.append(i.get_resname())
                non_aa_residues.append(i.get_resname())
        return aa_residues, seq, non_aa_residues
    
    def compute_distance_matrix(self, chain_1, chain_2):
        """
        compute distance matrix of two chains
        """
        dist_matrix = np.zeros((len(chain_1), len(chain_2)), np.float)
        for row, residue_1 in enumerate(chain_1):
            for col, residue_2 in enumerate(chain_2):
                dist_matrix[row, col] = self.compute_res_res_distance(
                    residue_1, residue_2)
        return dist_matrix    

    def compute_res_res_distance(self, residue_1, residue_2):
        """
        compute distance of two residue's alpha-carbon's coordinates
        """
        GLY = 'GLY'
        res_1_name = residue_1.get_resname()
        res_2_name = residue_2.get_resname()
        diff_vector = 0.0
        try:
            if res_1_name == GLY and res_2_name != GLY:
                diff_vector = residue_1["CA"].coord - residue_2["CB"].coord
            elif res_1_name != GLY and res_2_name == GLY:
                diff_vector = residue_1["CB"].coord - residue_2["CA"].coord
            elif res_1_name == GLY and res_2_name == GLY:
                diff_vector = residue_1["CA"].coord - residue_2["CA"].coord
            else:
                diff_vector = residue_1["CB"].coord - residue_2["CB"].coord
        except Exception as e:
            print("Can not resolve distance: ", res_1_name, res_2_name)
            # traceback.print_exc()
            raise

        return np.sqrt(np.sum(diff_vector * diff_vector))

    def get_contact_map(self, pdb_code, chain_id):
        print("computing contact-map for {}:{} ... ...".format(pdb_code, chain_id))
        pdb_filename = CONFIGS.PDB_DIR + pdb_code + CONFIGS.DOT_CIF
        is_defected = False
        # reading whole structure
        structure = self.parser.get_structure(pdb_code, pdb_filename)
        models = list(structure.get_models())
        chains = list(models[0].get_chains())
        # for each chain
        for chain in chains:
            if chain.id == chain_id:
                all_residues = list(chain.get_residues())
                aa_residues, seq, _ = self.filter_aa_residues(all_residues)
                n_aa_residues = len(aa_residues)
                dist_matrix = np.zeros((n_aa_residues, n_aa_residues), np.float)
                contact_map = np.zeros((n_aa_residues, n_aa_residues), np.float)
                try:
                    # computing distance matrix
                    dist_matrix = self.compute_distance_matrix(aa_residues, aa_residues)
                except Exception as e:
                    is_defected = True
                    break
                # computing comtact map on threshhold
                contact_map = np.where(dist_matrix < self.threshhold, 1, 0)
        
        DataUtils.save_contact_map(contact_map, pdb_with_chain)
        # DataViz.plot_images([contact_map], pdb_with_chain, cols=1)
        return is_defected, contact_map, dist_matrix

    def convert_cif_to_fasta(self, pdb_code, chain_id):
        """
        Bio.SeqIO describes how to convert in the following page in "Conversion" section
            https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
        The following link has an example for converting "cif to fasta" format.
            http://sequenceconversion.bugaco.com/converter/biology/sequences/cif-atom_to_fasta.php
        """
        print("computing fasta for {}:{} ... ...".format(pdb_code, chain_id))
        protein_cif = CONFIGS.PDB_DIR + pdb_code + CONFIGS.DOT_CIF
        records = SeqIO.parse(protein_cif, CONFIGS.CIF_ATOM)
        for record in records:
            if record.id[-1] == chain_id:
                protein_fasta = CONFIGS.FASTA_DIR + pdb_code + chain_id + CONFIGS.DOT_FASTA
                SeqIO.write(record, protein_fasta, CONFIGS.FASTA)

    def get_blast_result(self, pdb_code, chain_id):
        """
        From following link, section: "Running BLAST over the Internet" 
            http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc91
        """
        print("computing blastp for {} over internet".format(pdb_code))
        protein_fasta = CONFIGS.FASTA_DIR + pdb_code + chain_id + CONFIGS.DOT_FASTA
        record = SeqIO.read(protein_fasta, format=CONFIGS.FASTA) # read is for when you have exactly one object, and parse is an iterator for when you can have lots of objects 
        result_handle = NCBIWWW.qblast("blastp", "nt", record.seq)
        result_file = CONFIGS.BLAST_DIR + pdb_code + chain_id + CONFIGS.DOT_XML
        with open(result_file, "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

    def psi_blast(self, pdb_code, chain_id):
        """
        This blast run the query sequence 3 iterations against a db using psiblast program,
        and save the output file in psiblast directory as xml format. 
        """
        print("computing PSSM for {}:{} using psi-blast ... ...".format(pdb_code, chain_id))
        pdb_with_chain_id = pdb_code + chain_id
        query = CONFIGS.FASTA_DIR + pdb_with_chain_id + CONFIGS.DOT_FASTA
        # out_xml = CONFIGS.PSIBLAST_DIR + pdb_with_chain_id + CONFIGS.DOT_XML
        # out_pssm = CONFIGS.PSIBLAST_DIR + "pssm.txt"
        out_ascii_pssm = CONFIGS.PSIBLAST_DIR + pdb_with_chain_id + ".pssm"
        E_VALUE_TRESH = "10"
        
        cline = NcbipsiblastCommandline(cmd=self.psiblast_exe, db=self.blast_db, query=query,\
                                        evalue = E_VALUE_TRESH, outfmt = 5, num_iterations=3,\
                                        save_pssm_after_last_round=True, out_ascii_pssm=out_ascii_pssm)
                                        # out = out_xml, out_pssm=out_pssm, 
        cline()
        
    def run_scratch_suite(self, pdb_code, chain_id):
        """
        This function runs scratch-suite for each pdbid-chainid.
        SCRATCH currently generates four output files:
        .ss, .ss8, acc, acc20
        """
        print("running SCRATCH-suite for {}:{} ... ...".format(pdb_code, chain_id))
        input_fasta_file = "../fastas/" + pdb_code + chain_id + CONFIGS.DOT_FASTA
        command = "cd " + CONFIGS.SCRATCH_OUT + "; " + \
            CONFIGS.SCRATCH_EXE_PATH + " " + input_fasta_file + " " + pdb_code + chain_id + " " + str(4) +";"
        proc = subprocess.Popen(command, shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        proc.communicate() # comment me if you want to run the process in background


pdb = PDB()
pdb.run_scratch_suite("6y2d", "A")
# file_content = open(CONFIGS.ALL_PDB_IDS, "r")
# good_pdbs = []
# bad_pdbs = []
# for i, line in enumerate(file_content):
#     print("{}th protein:".format(i+1))
#     pdb_code, chain_id = pdb.read_pdb_ids(line)
#     pdb_with_chain = pdb_code + chain_id
#     # print(pdb_code, chain_id)
#     pdb.download(pdb_code)
#     is_defected, contact_map, dist_matrix = pdb.get_contact_map(pdb_code, chain_id)
#     pdb.convert_cif_to_fasta(pdb_code, chain_id)
#     pdb.psi_blast(pdb_code, chain_id)
#     pdb.run_scratch_suite(pdb_code, chain_id)
#     if not is_defected:
#         good_pdbs.append(pdb_with_chain)
#     else:
#         bad_pdbs.append(pdb_with_chain)
#     print()
# # save good_pdbs, and bad_pdbs in file
# DataUtils.save_itemlist(bad_pdbs, CONFIGS.BAD_PDB_IDS)
# DataUtils.save_itemlist(good_pdbs, CONFIGS.GOOD_PDB_IDS)
