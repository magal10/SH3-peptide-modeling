import subprocess as sp
import os
from Bio.PDB import PDBParser
import shutil

RMSD_TH = 2.0
BEST_MATCH_TH = 40
ALIGN_SCRIPT = "/cs/staff/dina/scripts/alignMP.pl"
CHAIN_SCRIPT = "/cs/staff/dina/scripts/extractChains.pl"
TRANS_SCRIPT = "/cs/staff/dina/utils/pdb_trans"


def align(directory, dataname):
    ROOT_DIR = os.path.abspath(os.curdir)
    model = "data/modelSH3.pdb"
    res = open("Newresults" + dataname + ".txt", "w")
    filelist=os.listdir(directory)
    for filename in filelist:  # goes through PDB files
        f = os.path.join(directory, filename)
        if os.path.isfile(f) and f.endswith(".pdb"):
            name = filename.split('.')[0]
            os.mkdir(directory + "/" + name + "dir")
            shutil.copy("params.txt" , directory + "/" + name + "dir")

            os.chdir(directory + "/" + name + "dir")

            sp.run([CHAIN_SCRIPT, os.path.join("../", filename)])

            # best
            best_RMSD = 0
            best_match_size = 0
            best_chain = ""
            best_transformation = ""

            # current
            current_RMSD = 0
            current_match_size = 0
            current_transformation = ""
            chainfiles = os.listdir(os.getcwd())
            for file in chainfiles:
                chain = file.split(".")[1]
                sp.run([ALIGN_SCRIPT, "../../../" + model, file, name + chain + ".pdb"])
                if os.path.exists("2_sol.res"):
                    for line in open("2_sol.res"):
                        if line.startswith('Match List'):
                            current_match_size = int(line.split()[-1])

                        if line.startswith('RMSD : '):
                            current_RMSD = float(line.split()[-1])

                        if line.startswith('Trans : '):
                            current_transformation = line.split(":")[1].strip()


                    if current_match_size > BEST_MATCH_TH and current_RMSD < RMSD_TH:
                        if current_match_size > best_match_size:
                            best_match_size = current_match_size
                            best_RMSD = current_RMSD
                            best_chain = chain
                            best_transformation = current_transformation


                        elif current_match_size == best_match_size:
                            if current_RMSD > best_RMSD:
                                best_match_size = current_match_size
                                best_RMSD = current_RMSD
                                best_chain = chain
                                best_transformation = current_transformation

            os.chdir(ROOT_DIR)
            res.write(filename + "\t" + best_chain + "\t" + best_transformation + "\n")

