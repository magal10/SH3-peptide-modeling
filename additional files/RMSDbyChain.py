import subprocess as sp
import os

RMSD_SCRIPT = "/cs/staff/dina/utils/rmsd"
TRANS_SCRIPT = "/cs/staff/dina/utils/pdb_trans"
EXTRACT_CHAIN_SCRIPT = "/cs/staff/dina/utils/getChain.Linux"


def rmsdByChain(directory, predictor, predDir,domainChain, peptideChain ):
    rmsdFile = open("RMSD" + predictor + ".txt", "w")
    rmsdBoxPlot = open("RMSDforBoxplot" + predictor + ".txt", "w")

    filelist = os.listdir(directory)
    for filename in filelist:  # goes through PDB files
        f = os.path.join(directory, filename)
        if os.path.isfile(f) and f.endswith(".pdb"):
            ref = predDir + "/" + filename
            model = directory + "/" + filename
            modelSh3 = directory + "/test_sh3/sh3" + filename
            modelTransformed = directory + "/testTransformed/" + filename

            modelTransformedA = directory + "/testTransformed/A" + filename
            modelTransformedB = directory + "/testTransformed/B" + filename

            refSh3 = predDir + "/sh3/" + filename
            refPeptide = predDir + "/peptide/" + filename

            # align reference SH3 to the predicted SH3
            sp.run(f"{EXTRACT_CHAIN_SCRIPT} A {model} > {modelSh3}", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            sp.run(f"{EXTRACT_CHAIN_SCRIPT} {domainChain} {ref} > {refSh3}", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            transformation = str(sp.run(f"{RMSD_SCRIPT} {refSh3} {modelSh3} -t | head -n1",
                                        shell=True, capture_output=True).stdout.strip().decode("utf-8"))

            # transform the reference structure by the transformation calculated above
            sp.run(f"{TRANS_SCRIPT} {transformation} < {model} > {modelTransformed}", shell=True, stdout=sp.DEVNULL,
                   stderr=sp.DEVNULL)

            # divide PDBs by chains
            sp.run(f"{EXTRACT_CHAIN_SCRIPT} {peptideChain} {ref} > {refPeptide}", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            sp.run(f"{EXTRACT_CHAIN_SCRIPT} A {modelTransformed} > {modelTransformedA}", shell=True, stdout=sp.DEVNULL,
                   stderr=sp.DEVNULL)
            sp.run(f"{EXTRACT_CHAIN_SCRIPT} B {modelTransformed}  > {modelTransformedB}", shell=True, stdout=sp.DEVNULL,
                   stderr=sp.DEVNULL)

            # calculate RMSDs
            rmsdSh3 = float(sp.run(f"{RMSD_SCRIPT} {refSh3}  {modelTransformedA}| tail -n1 ", shell=True,
                                   capture_output=True).stdout.strip())
            rmsdPeptide = float(sp.run(f"{RMSD_SCRIPT} {refPeptide} {modelTransformedB}| tail -n1 ", shell=True,
                                       capture_output=True).stdout.strip())

            # write to files
            rmsdFile.write(filename + "\t" + str (rmsdSh3) + "\t" + str (rmsdPeptide) + "\n")
            rmsdBoxPlot.write(str (rmsdSh3) + "\tdomain\n" + str (rmsdPeptide) + "\tpeptide\n")

