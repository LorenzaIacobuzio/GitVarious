#!/usr/bin/env python

import os
import sys
import re
import optparse
import argparse
import glob
import shutil
import itertools
import stat
import mmap
import fnmatch
import tempfile
from itertools import izip_longest

usage = "%RunJobs.py Binary [options]"
parser = optparse.OptionParser(usage)
arg = argparse.ArgumentParser(usage)

# Arg parser

parser.add_option("-i", "--info", help="Script to submit jobs for analysis; Binary is analyzer name, List contains the reco files to be analyzed")
parser.add_option("-l", "--launch", dest="launch", action="store_true", default=False, help="if enabled, launches jobs")
parser.add_option("-a", "--add", dest="add", action="store_true", default=False, help="if enabled, hadds output histograms of each job")
parser.add_option("-q", "--queue", dest="queue", default="2nd", help="type of queue; if not specified, default is 2nd")
parser.add_option("-p","--path", dest="path", default="/afs/cern.ch/na62/offline/lists/Data/Reco/2016/", help="path of files to be processed; if not specified, default is /afs/cern.ch/na62/offline/lists/Data/Reco/2016/")
parser.add_option("-o","--output", dest="output", default="/afs/cern.ch/work/l/liacobuz/private/", help="path of output files; if not specified, default is /afs/cern.ch/work/l/liacobuz/private/")

(options, args) = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

arg.add_argument('args', nargs='*')

binary = os.path.abspath(sys.argv[1])
launch = options.launch
add = options.add
queue = options.queue
path = options.path
output = options.output
listpath = "/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/Lists/Lists/"
jobpath = "/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/Lists/Job/"

# Function for creating single job script

def JobScript(binary, listname, output):
    with open(jobpath + "singleJob" + listname + ".sh", "w") as f:
        f.write("#!/usr/bin/env bash \ncd /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/ \nsource /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh \n" + binary + " -l " + listpath + listname + " -o " + output + listname + ".root --ignore -n 100")
    os.system("chmod u+x " + jobpath + "singleJob" + listname + ".sh")

# Function for splitting reco files into sublists

def JobManager(binary, path, queue):
    ListOfLists = os.listdir(path)
    NLists = len(ListOfLists)
    NFilesPerList = 150
    print("[RunJobs] Removing existing sublists from " + listpath)
    for file in os.listdir(listpath):
        if "list" in file:
            os.remove(os.path.join(listpath, file))
    print("[RunJobs] Creating new list of all reco files in " + path)
    with open(listpath + "/List", "w") as f:
        for file in ListOfLists:
            if "HNLVtx" in file and "r1666" not in file and "EOS" in file and "v0.11.1" in file:
                with open(str(path + file), "r") as d:
                    for line in d:
                        f.write(line)
    print("[RunJobs] Splitting into sublists in " + listpath)
    with open(listpath + "List", "r") as f:
        os.system("split -l " + str(NFilesPerList) + " " + listpath + "List " + listpath + "Sublist-")

# Execution

counter = 0

if launch:
    print(str("[RunJobs] About to execute binary " + binary + ". Output file will be in " + output))
    JobManager(binary, path, queue)
    print(str("[RunJobs] Deleting job scripts in " + jobpath))
    for file in os.listdir(jobpath):
        os.remove(os.path.join(jobpath, file))
    for f in os.listdir(listpath):
        if "Sublist" in f and counter < 3:
            counter += 1
            JobScript(binary, f, output)
            os.system("bsub -q " + queue + " -oo " + output + "log-" + f + " " + jobpath + "singleJob" + f + ".sh")

if add:
    print(str("[RunJobs] About to hadd all output histograms. Output will be in " + output))
    subprocess.call(str("hadd " + output + "Result.root " + output + "*.root"))
