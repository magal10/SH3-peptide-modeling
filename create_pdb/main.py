import csv
import subprocess
import sys

if __name__ == '__main__':

    peptideDic = {}
    sh3Dic ={}
    with open(sys.argv[1],mode="r") as fp:
        readerP = csv.reader(fp)
        peptideDic = {rows[0]:rows[1:4] for rows in readerP}

    with open(sys.argv[2],mode="r") as fs:
        readerS = csv.reader(fs)
        sh3Dic = {rows[0]:rows[1:4] for rows in readerS}

    max_len_p = -1
    max_len_s = -1


    for key in peptideDic:

        script_get_frag = "/cs/staff/dina/utils/get_frag_chain.Linux"
        arg_original_pdb = "/cs/usr/magal.gafni/PycharmProjects/create_pdb/transformedPDBs/best" + key
        arg_chain_p = peptideDic[key][0]
        arg_start_p = peptideDic[key][1]
        arg_end_p = peptideDic[key][2]
        outfile_p = "/cs/usr/magal.gafni/PycharmProjects/create_pdb/peptide_pdbs/P"+key

        subprocess.run(f"{script_get_frag} {arg_original_pdb} {arg_chain_p} {arg_start_p} {arg_end_p} > {outfile_p}",shell=True)
        subprocess.run(["/cs/staff/dina/scripts/namechain.pl","/cs/usr/magal.gafni/PycharmProjects/create_pdb/peptide_pdbs/P"+key, "B" ])


        arg_chain_s = sh3Dic[key][0]
        arg_start_s = sh3Dic[key][1]
        arg_end_s = sh3Dic[key][2]
        outfile_s = "/cs/usr/magal.gafni/PycharmProjects/create_pdb/sh3_pdbs/S"+key

        subprocess.run(f"{script_get_frag} {arg_original_pdb} {arg_chain_s} {arg_start_s} {arg_end_s} > {outfile_s}",shell=True)
        subprocess.run(["/cs/staff/dina/scripts/namechain.pl","/cs/usr/magal.gafni/PycharmProjects/create_pdb/sh3_pdbs/S"+key, "A" ])

        script_join = "/cs/staff/dina/scripts/joinModelsPDB.pl"
        pdb1 = "sh3_pdbs/S"+key
        pdb2 = "peptide_pdbs/P"+key
        outfile = "ready_pdbs/" +key

        subprocess.run(f"cat {pdb1} {pdb2} > {outfile}",shell = True)

        if int(peptideDic[key][2]) - int(peptideDic[key][1]) > max_len_p:
            max_len_p = int(peptideDic[key][2]) - int(peptideDic[key][1])

        if int(sh3Dic[key][2]) - int(sh3Dic[key][1]) > max_len_s:
            max_len_s = int(sh3Dic[key][2]) - int(sh3Dic[key][1])

    print( "max len pep: "+str(max_len_p))
    print("max len sh3: " + str(max_len_s))



