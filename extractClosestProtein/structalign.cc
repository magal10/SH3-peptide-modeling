#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <chrono>
#include <iostream>
#include "Triangle.h"

int main(int argc , char* argv[]){

    if(argc !=4) {
        std::cerr << "Usage: "<<argv[0]<< " dist_threshold target.pdb model.pdb" << std::endl;
        exit(1);
    }

    //********Parameters********************
    float m_fDistThr = atof(argv[1]); // distance threshold on atoms in correspondence


    // read the two files into Molecule
    Molecule<Atom> molModel, molTarget, molTemp;

    std::ifstream fileModel(argv[3]);
    std::ifstream fileTarget(argv[2]);
    std::ifstream fileTemp(argv[3]);

    if(!fileModel) {
        std::cout<< "File " << argv[3] << "does not exist." << std::endl;
        return 0;
    }
    if(!fileTarget) {
        std::cout << "File " << argv[2] << " does not exist." << std::endl;
        return 0;
    }

    molTemp.readPDBfile(fileTemp);
    if(molTemp[0].isRNABackbone()){
        molModel.readPDBfile(fileModel, PDB::PSelector());
        molTarget.readPDBfile(fileTarget, PDB::PSelector());
    }else{
        molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
        molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
    }


    // next we insert the target molecule into hash
    // this will help us to find atoms that are close faster
    GeomHash <Vector3,int> gHash(3,m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
    for(unsigned int i=0; i<molTarget.size(); i++) {
        gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
    }

    // now we try random rotations and choose the best alignment from random rotations
    unsigned int iMaxSize=0;
    RigidTrans3 rtransBest;
    float rmsd;
    for(unsigned int i=0; i< molModel.size()-2; i++){
        Vector3 model_atom1 = molModel[i].position();
        Vector3 model_atom2 = molModel[i+1].position();
        Vector3 model_atom3 = molModel[i+2].position();
        Triangle trModel = Triangle(model_atom1,model_atom2,model_atom3);
        for(unsigned int j=0; j< molTarget.size()-2; j++){
            Vector3 target_atom1 = molTarget[j].position();
            Vector3 target_atom2 = molTarget[j+1].position();
            Vector3 target_atom3 = molTarget[j+2].position();
            Triangle trTarget = Triangle(target_atom1,target_atom2,target_atom3);

            //find the transformation for the triangles
            RigidTrans3 curTrans = trTarget | trModel;

            // match is a class that stores the correspondence list, eg.
            // pairs of atoms, one from each molecule, that are matching
            Match match;

            // apply rotation on each atom in the model molecule and
            // add the pairs of atoms (one from target and one from model)
            // that are close enough to the match list
            for(unsigned int t=0; t< molModel.size(); t++) {
                Vector3 transformed_atom = curTrans*molModel[t].position(); // rotate

                // find close target molecule atoms using the hash
                HashResult<int> result;
                gHash.query(transformed_atom, m_fDistThr, result); // key is mol atom coordinate

                // check if the atoms in the result are inside the distance threshold
                // the hash is a cube shape, there can be atoms further that the threshold
                for(auto x = result.begin(); x != result.end(); x++) {
                    float dist = transformed_atom.dist(molTarget[*x].position());
                    if(dist <= m_fDistThr) {
                        float score = (1 / (1 + dist));
                        match.add( *x , t, score, score );
                    }
                }
                result.clear();
            }

            //calculates transformation that is a little better than "rotation"
            match.calculateBestFit(molTarget, molModel);

            if(iMaxSize < match.size() ){
                iMaxSize = match.size();
                rtransBest=match.rigidTrans();
                rmsd = match.rmsd();
            }
        }

    }

    std::cout << iMaxSize << " " << rmsd << " " << rtransBest.translation() << " " <<
              rtransBest.rotationAngles() << std::endl;

    std::ofstream output_file;
    molTemp *= rtransBest;
    output_file.open("transformed.pdb");
    output_file << molTemp;
    output_file.close();

}
