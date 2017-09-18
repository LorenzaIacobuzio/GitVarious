#!/usr/bin/env python

import os
import sys
import optparse
import argparse
import shutil
import stat
import time
import re
import glob

usage = "%CheckJobs.py Binary List [options]"
parser = optparse.OptionParser(usage)
arg = argparse.ArgumentParser(usage)

# Arg parser

parser.add_option("-i", "--info", help="List is a text file with the run numbers to process on CASTOR/EOS; either option -l or -a must be specified")
parser.add_option("-l", "--launch", dest="launch", action="store_true", default=False, help="Submit jobs")
parser.add_option("-a", "--add", dest="add", action="store_true", default=False, help="Hadd output histograms")
parser.add_option("-q", "--queue", dest="queue", default="8nh", help="type of queue; if not specified, default is 8nh")
parser.add_option("-j","--job", dest="job", default="job", help="job name; if not specified, default is job")
parser.add_option("-n","--number", dest="number", default="30", help="number of files per job; if not specified, default is 30")
parser.add_option("-s","--storage", dest="storage", default="EOS", help="storage of files to process; if not specified, default is EOS")
parser.add_option("-p","--path", dest="path", default="/afs/cern.ch/na62/offline/lists/Data/Reco/2016/", help="path of files (on CASTOR/EOS) to be processed; if not specified, default is /afs/cern.ch/na62/offline/lists/Data/Reco/2016/")
parser.add_option("-o","--output", dest="output", default="./MyAnalyzerOutput", help="path of output files; if not specified, default is ./MyAnalyzerOutput")

(options, args) = parser.parse_args()

if len(sys.argv) < 4:
    parser.print_help()
    sys.exit(1)

arg.add_argument('args', nargs='*')

launch = options.launch
add = options.add
binary = os.path.abspath(sys.argv[1])
binary_name = os.path.split(sys.argv[1])[1]
list_of_runs = os.path.abspath(sys.argv[2])
queue = options.queue
job = options.job
number = options.number
storage = options.storage
path = options.path
output = "/".join(os.path.abspath(options.output).split('/')[:-1])
name = os.path.basename(os.path.normpath(options.output))

list_of_files = os.listdir(path)    
files_in_dir = len(list_of_files)

# Function for job handling

def job_manager(binary, binary_name, number, job, queue):
    counter = 0;
    prefix = "xfarmer"
    os.system(str("split -a 3 -l " + number + " " + output + "/FullListOfRuns.txt " + prefix))
    current_dir = os.getcwd()
    farmer_in_run = os.listdir(current_dir)
    nj = 0
    for item in farmer_in_run:
        if prefix not in item:
            continue
        elif prefix in item:
            if counter < 300:
                item.rstrip()
                directory = str("x" + item)
                os.makedirs(directory)
                os.symlink(str(binary), str(current_dir + "/" + directory + "/" + binary_name))
                shutil.move(str(item), str(directory + "/list"))
                os.chdir(str(current_dir + "/" + directory))
                njs = nj
                if len(str(njs)) < 2:
                    njsf = "{0:0>2}".format(njs)
                else:
                    njsf = str(njs)
                njfull = job + njsf
                with open("job", "w") as f:
                    f.write(str("#!/bin/bash\ncd " + current_dir + "/" + directory + "\nsource /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh\n" + current_dir + "/" + directory + "/" + binary_name + " -l list -o outFile.root\n"))
                st = os.stat("job")
                os.chmod("job", st.st_mode | stat.S_IEXEC)
                os.system(str("bsub -q " + queue + " -o " + current_dir + "/" + directory + "/log -e " + current_dir + "/" + directory + "/err -J " + njfull + " " + current_dir + "/" + directory + "/job"))
                os.chdir(current_dir)
                nj += 1
                counter += 1
            elif counter == 300:
                with open(str(output + "/Stop.txt"), "w") as f:
                    with open(str(output + "/FullListOfRuns.txt"), "r") as g:
                        stopline = g.readlines()[300*int(number)].rstrip()
                        f.write(str("Stopped after processing file " + stopline))
                        print(str("Stopped after processing file " + stopline + ". No more than 300 jobs allowed"))                                        
                        counter += 1
            elif counter > 300:
                if re.search("xfarmer", item):
                    print("Cleaning up...")
                    os.remove(os.path.join(os.getcwd(), item))
    for item in os.listdir(os.getcwd()):
        if re.search("^xfarmer", item):
            os.remove(os.path.join(os.getcwd(), item))
    sys.exit(1)
    return

# Create joined list of files

if launch:
    os.chdir(output)
    with open(list_of_runs, "r") as f:
        with open(str(output + "/FullListOfRuns.txt"), "w") as g:
            for line in f:
                useless_files = 0
                for file in list_of_files:
                    if line.rstrip() in file and "Electron" not in file and "3Tracks" not in file and storage.upper() in file:
                        file_name = str(path + file)
                    elif line.rstrip() not in file:
                        useless_files += 1
                if files_in_dir == useless_files:
                    print(str("[Warning]: Run" + line.rstrip() + " not found"))
                else:
                    with open(file_name, "r") as h:
                        if os.stat(file_name).st_size == 0:
                            print(str("[Warning]: Run" + line.rstrip() + " is an empty file"))
                        else:
                            s = h.read()
                            g.write(s)

# Submit jobs

    print(str("[CheckJobs] About to execute binary " + binary_name + ". Output file is " + output + "/" + name))                                    
    if not os.path.exists(os.path.abspath(name)):
        os.makedirs(str(output + "/" + name))
    elif os.path.exists(os.path.abspath(name)):
        print("Directory for runs already exists! Remove it before executing this script")
        sys.exit(1)
    os.chdir(str(output + "/" + name))
    job_manager(binary, binary_name, number, job, queue)
    os.chdir(str(output + "/" + name))

# Hadd histograms

if add:
    print(str("[CheckJobs] About to hadd all histograms. Output file is in " + output + "/" + name))
    os.chdir(str(output + "/" + name))
    if os.path.exists(str("Result" + binary_name + ".root")):
        os.remove(str("Result" + binary_name + ".root"))
    histos_to_add = " ".join(glob.glob('./*/outFile.root'))
    os.system(str("hadd Result" + binary_name + ".root " + histos_to_add))
