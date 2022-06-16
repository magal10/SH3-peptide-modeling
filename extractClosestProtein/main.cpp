#include <iostream>
#include "Vector3.h"
#include "Atom.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <map>
#include <string>
#include <boost/algorithm/string.hpp>

#define PEPTIDE "peptide"
#define DOMAIN "domain"

std::map<std::string, char> readToMap(char* path,std::map<std::string, char> sh3ChainMap){
    std::ifstream sh3Chains(path);
    if(!sh3Chains){
        std::cerr << "no such file 3" << std::endl;
    }
    std::string line;
    while (!sh3Chains.eof())
    {
        getline(sh3Chains, line);
        boost::trim(line); // remove all spaces
        if (line.length() ==0 )
            continue;
        // skip comments
        if (line[0] == '#' || line[0] == '\0')
            continue;

        std::vector<std::string> split_results;
        boost::split(split_results, line, boost::is_any_of("\t ,"),
                     boost::token_compress_on);
        if(split_results.size() !=2){
            continue;
        }
        std::string key = split_results[0];
        char val = split_results[1][0];
        sh3ChainMap[key] = val;
    }
    sh3Chains.close();
    return sh3ChainMap;
}

int main(int argc , char* argv[]){

    unsigned int peptideTooLong = 0 ;
    unsigned int noSuchFile = 0;
    int no_peptide = 0;
    float m_fDistThr = 5;
    int loop_counter =0;
    std::string proteinTarget(argv[3]);
    Molecule<Atom> reference;
    std::ifstream fileReference(argv[1]);
    if(!fileReference){
        std::cerr << "no such file 2"  << std::endl;
    }
    reference.readPDBfile(fileReference, PDB::CAlphaSelector());

    std::map<std::string, char> sh3ChainMap;
    sh3ChainMap = readToMap(argv[2],sh3ChainMap);

    std::fstream csvOut;
    csvOut.open(proteinTarget+"_chain_res.csv", std::ios::out| std::ios::app);


    for(auto iter = sh3ChainMap.begin() ; iter!= sh3ChainMap.end(); ++iter) {
        loop_counter++;
        std::string pdbName = iter->first;

        std::cout << "start iterate pdbs: " << pdbName << std::endl;
        Molecule<Atom> molModel;
        std::ifstream fileModel("transformedPDBs/best"+pdbName);
        if(!fileModel){
            std::cerr << "no such file " << pdbName<< std::endl;
            noSuchFile++;
            continue;
        }
        molModel.readAllPDBfile(fileModel,PDB::CAlphaSelector());

        GeomHash <Vector3,int> gHash(3,m_fDistThr);

        for(unsigned int i=0; i<molModel.size(); i++) {
            if(proteinTarget == PEPTIDE){
                if(molModel[i].chainId()!= sh3ChainMap[pdbName]){ //ONLY FOR PEPTIDE
                    gHash.insert(molModel[i].position(), i); // coordinate is the key to the hash, we store atom index
                }
            }
            else{
                if(molModel[i].chainId() == sh3ChainMap[pdbName]){ //ONLY FOR SH3
                    gHash.insert(molModel[i].position(), i); // coordinate is the key to the hash, we store atom index
                }
            }

        }
        HashResult<int> result;
        Match match;
        for (int i = 0; i < reference.size(); ++i) {
            gHash.query(reference[i].position(), m_fDistThr, result);

            for(auto x = result.begin(); x != result.end(); x++) {
                float dist = reference[i].position().dist(molModel[*x].position());
                if(dist <= m_fDistThr) {
                    float score = (1 / (1 + dist));
                    match.add( *x , i, score, score );
                }
            }
            result.clear();

        }
        std::map<char, std::vector<unsigned int>> potential;
        std::cout << pdbName << " match size: " << match.size() << std::endl;

        for (int i = 0; i <match.size() ; ++i) {
            unsigned int molModelIndex = match[i].model;
            unsigned int resIndex = molModel[molModelIndex].residueIndex();
            char chainId = molModel[molModelIndex].chainId();
            if(potential.find(chainId) == potential.end()){
                potential[chainId] = {resIndex,resIndex,1};
            }
            else{
                potential[chainId][2]++;
                if(resIndex < potential[chainId][0]){
                    potential[chainId][0] = resIndex;
                }
                if(resIndex > potential[chainId][1]){
                    potential[chainId][1] = resIndex;
                }
            }
        }
        if(potential.empty()){
            no_peptide++;
            continue;
        }
        char maxChainID;
        unsigned int maxRes = 0;
        for(auto x = potential.begin(); x != potential.end(); x++){
            if(x->second[2] > maxRes){
                maxChainID= x->first;
                maxRes = x->second[2];
            }
        }


        if(potential[maxChainID][1] - potential[maxChainID][0] > 20 && proteinTarget==PEPTIDE){ //ONLY FOR PEPTIDE
            peptideTooLong++;
            continue;
        }
        csvOut << pdbName << ", " << maxChainID << ", " << potential[maxChainID][0] << ", " //ONLY FOR PEPTIDE
                << potential[maxChainID][1] << "\n";

    }
        std::cout << "file without peptide : " << no_peptide << std::endl;
        std::cout << "file with peptide too long : " << peptideTooLong << std::endl;
        std::cout << "loop counter : " << loop_counter << std::endl;
        std::cout << "mapSize : " << sh3ChainMap.size() << std::endl;
        csvOut.close();

}

