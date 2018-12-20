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

parser.add_option("-i", "--info", help = "Script to submit jobs for analysis; Binary is the analyzer name")
parser.add_option("-q", "--queue", dest = "queue", default = "1nd", help = "type of queue; default is 1nd")
parser.add_option("-p", "--path", dest = "path", default = "/afs/cern.ch/na62/offline/lists/Data/Reco/2016A-v1.0.0/", help = "path of files to be processed; default is /afs/cern.ch/na62/offline/lists/Data/Reco/2016A-v1.0.0/")
parser.add_option("-o", "--output", dest = "output", default = "/afs/cern.ch/work/l/liacobuz/private/", help = "path of output files; default is /afs/cern.ch/work/l/liacobuz/private/")
parser.add_option("-n", "--number", dest = "number", default = "50", help = "number of files per sublist; default is 50")

(options, args) = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

arg.add_argument('args', nargs='*')

binary = os.path.abspath(sys.argv[1])
queue = options.queue
path = options.path
output = options.output
number = options.number
listpath = output + "Lists/Lists/"
jobpath = output + "Lists/Job/"
binarypath = (binary.rpartition("/")[0]).rpartition("/")[0] + "/"

# Function for creating single job script

def JobScript(binary, listname, output):
    with open(jobpath + "singleJob" + listname + ".sh", "w") as f: # creating job script with sublist name
        f.write("#!/usr/bin/env bash \ncd " + binarypath + " \nsource " + binarypath + "scripts/env.sh \n" + binary + " -l " + listpath + listname + " -o " + output + listname + ".root --ignore") # script executes analyzer on sublist
    os.system("chmod u+x " + jobpath + "singleJob" + listname + ".sh") # permissions for script

# Function for splitting reco files into sublists

def JobManager(binary, path, queue):
    ListOfLists = os.listdir(path) # listing all reco lists
    NLists = len(ListOfLists)

    if os.path.exists(listpath): # if sublist directory exists, empty it
        print("[RunJobs] Removing existing sublists from " + listpath)
        for file in os.listdir(listpath):
            if "list" in file:
                os.remove(os.path.join(listpath, file))
    else: # otherwise, create it
        os.makedirs(listpath)

    print("[RunJobs] Creating new list of all reco files in " + path) # creating list of all reco files in reco lists

    with open(listpath + "/List", "w") as f:
        for file in ListOfLists:
            if "HNL" in file and "EOS" in file and "v1.0.0" in file:
                with open(str(path + file), "r") as d:
                    for line in d:
                        f.write(line)

    print("[RunJobs] Splitting into sublists in " + listpath) # splitting list into sublists

    with open(listpath + "List", "r") as f:
        os.system("split -a 5 -l " + str(number) + " " + listpath + "List " + listpath + "Sublist-")

# Main

def main():
    #counter = 0

    print(str("[RunJobs] About to execute binary " + binary + ". Output file will be in " + output))

    JobManager(binary, path, queue)

    if os.path.exists(jobpath): # if job script directory exists, empty it
        print(str("[RunJobs] Deleting job scripts in " + jobpath))
        for file in os.listdir(jobpath):
            os.remove(os.path.join(jobpath, file))
    else: # otherwise, create it
        os.makedirs(jobpath)

    for f in os.listdir(listpath):
        #if "Sublist" in f and counter < 1:
        if "Sublist" in f:
            #counter += 1
            JobScript(binary, f, output) # creating job script for specific sublist
            os.system("bsub -q " + queue + " -oo " + output + "log-" + f + " " + jobpath + "singleJob" + f + ".sh") # submitting specific job

# Execution

if __name__ == "__main__":
    main()
