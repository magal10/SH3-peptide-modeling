cmake_minimum_required(VERSION 3.20)
project(peptide)

set(CMAKE_CXX_STANDARD 17)

add_executable(peptide main.cpp Atom.h Atom.cc GeomHash.h HashLink.h HashResult.h macros.h Match.h Match.cc Matrix3.cc
        Matrix3.h Molecule.h numerics.cc numerics.h PDB.cc PDB.h RigidTrans3.cc RigidTrans3.h AminoAcid.h Rotation3.cc Rotation3.h
        Triangle.cc Triangle.h Vector3.cc Vector3.h)
