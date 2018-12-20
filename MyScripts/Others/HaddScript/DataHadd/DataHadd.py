#!/usr/bin/env python

# before, if copyFiles is not compiled:                                   
# C1: result of command root-config --cflags                         
# C2: result of command root-config --glibs                         
# g++ -o outfile C1 C2 copyFile.C         
# (on lxplus this means running g++ -o copyFiles -pthread -std=c++1y -m64 -I/cvmfs/sft.cern.ch/lcg/releases/LCG_86/ROOT/6.08.00/x86_64-slc6-gcc49-opt/include -L/cvmfs/sft.cern.ch/lcg/releases/LCG_86/ROOT/6.08.00/x86_64-slc6-gcc49-opt/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic copyFiles.C)   

# FOR DATA USE ONLY

import os
import sys
import argparse
import subprocess

usage = "%Hadd.py condorpath analyzer copypath"
parser = argparse.ArgumentParser(description = "Script for customized hadd. It calls copyFiles on each .root produced by Condor, creates a new .root with only the selected directory and hadds all of them in a final .root file")

parser.add_argument("condorpath", metavar = "condorpath", type = str, help = "leads to the condor outFiles directory")
parser.add_argument("analyzer", metavar = "analyzer", type = str, help = "analyzer directory name to be kept in the new .root file")
parser.add_argument("copypath", metavar = "copypath", type = str, help = "path to copyFiles.C")

if len(sys.argv) != 4:
    parser.print_help()
    sys.exit(1)

condorpath = os.path.abspath(sys.argv[1])
analyzer = sys.argv[2]
copypath = os.path.abspath(sys.argv[3])

def main():

    for subdir, dirs, files in os.walk(condorpath):
        for file in files:
            if "root" in file and "final" not in file and "new" not in file:
                os.system(copypath + "/copyFiles " + os.path.abspath(os.path.join(subdir, file)) + " " + os.path.abspath(os.path.join(subdir, "new" + file)) + " " + analyzer + "\n")
                hadd += " " + os.path.abspath(os.path.join(subdir, "new" + file))

    os.system(hadd)

    for subdir, dirs, files in os.walk(condorpath):
        for file in files:
            if "new" in file and "final" not in file:
                os.remove(os.path.abspath(os.path.join(subdir, file)))

if __name__ == "__main__":
    main()
