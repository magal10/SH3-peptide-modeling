import subprocess as sp
import os


RESULTS = "NewresultsNMR.txt"
TRANS_SCRIPT = "/cs/staff/dina/utils/pdb_trans"


def transform():
    with open(RESULTS) as res:
        for line in res:
            line = line.strip()
            line_arr = line.split("\t")
            if len(line_arr) == 3:

                pdb_name = "data/NMR/" + line_arr[0]
                print(pdb_name)
                chain = line_arr[1]
                transformation = line_arr[2]
                print(os.getcwd())
                output = "best"+line_arr[0]
                print(output)

                sp.run(f"{TRANS_SCRIPT} {transformation} < {pdb_name} > {output}", shell=True, stdout=sp.DEVNULL,
                       stderr=sp.DEVNULL)




