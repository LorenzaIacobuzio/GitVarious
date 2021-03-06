#!/usr/bin/env python

# before, if copyFiles is not compiled:                                   
# C1: result of command root-config --cflags                         
# C2: result of command root-config --glibs                         
# g++ -o outfile C1 C2 copyFiles.C         
# (on lxplus this means running g++ -o copyFiles -pthread -std=c++1y -m64 -I/cvmfs/sft.cern.ch/lcg/releases/LCG_86/ROOT/6.08.00/x86_64-slc6-gcc49-opt/include -L/cvmfs/sft.cern.ch/lcg/releases/LCG_86/ROOT/6.08.00/x86_64-slc6-gcc49-opt/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic copyFiles.C)

# FOR MC USE ONLY

import os
import sys
import argparse
import subprocess

usage = "%Hadd.py condorpath analyzer1 analyzer2 tree1 tree2 copypath"
parser = argparse.ArgumentParser(description = "Script for customized hadd. It calls copyFiles on each .root produced by Condor, creates a new .root with only the selected directory and hadds all of them in a final .root file")

parser.add_argument("condorpath", metavar = "condorpath", default = "/afs/cern.ch/work/l/liacobuz/private/condorOutput/HNL/outFiles/", type = str, help = "leads to the condor outFiles directory")
parser.add_argument("analyzer1", metavar = "analyzer1", type = str, help = "analyzer directory name to be kept in the new .root file")
parser.add_argument("analyzer2", metavar = "analyzer2", type = str, help = "analyzer directory name to be kept in the new .root file")
parser.add_argument("tree1", metavar = "tree1", type = str, help = "tree name to be kept in the new .root file")
parser.add_argument("tree2", metavar = "tree2", type = str, help = "tree name to be kept in the new .root file")
parser.add_argument("copypath", metavar = "copypath", type = str, help = "path to copyFiles.C")

if len(sys.argv) != 7:
    parser.print_help()
    sys.exit(1)

condorpath = os.path.abspath(sys.argv[1])
analyzer1 = sys.argv[2]
analyzer2 = sys.argv[3]
tree1 = sys.argv[4]
tree2 = sys.argv[5]
copypath = os.path.abspath(sys.argv[6])

def main():

    # Split files into subdirs (do only once at beginning)

    subprocess.check_call(['/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/HaddScript/SplitIntoSubfolders.sh', condorpath])

#    for subdir, dirs, files in os.walk(condorpath):
#        for file in files:
#            if "root" in file and "final" not in file and "new" not in file:
#                os.system(copypath + "/copyFiles " + os.path.abspath(os.path.join(subdir, file)) + " " + os.path.abspath(os.path.join(subdir, "new1" + file)) + " " + analyzer1 + " " + tree1 + "\n")
#                os.system(copypath + "/copyFiles " + os.path.abspath(os.path.join(subdir, file)) + " " + os.path.abspath(os.path.join(subdir, "new2" + file)) + " " + analyzer2 + " " + tree2 + "\n")
#
#    subprocess.call(['/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/HaddScript/MCHadd.sh', condorpath])

#    for subdir, dirs, files in os.walk(condorpath):
#        for file in files:
#            if "new" in file and "final" not in file and "Final" not in file:
#                os.remove(os.path.abspath(os.path.join(subdir, file)))

if __name__ == "__main__":
    main()
