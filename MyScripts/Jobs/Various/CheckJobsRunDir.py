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

usage = "%CheckJobs.py Binary List [options]"
parser = optparse.OptionParser(usage)
arg = argparse.ArgumentParser(usage)

# Arg parser

parser.add_option("-i", "--info", help="List is a text file with the run numbers to process on CASTOR/EOS")
parser.add_option("-l", "--launch", dest="launch", action="store_true", default=False, help="if enabled, launches jobs")
parser.add_option("-a", "--add", dest="add", action="store_true", default=False, help="if enabled, hadds output histograms of each job")
parser.add_option("-r", "--run", dest="run", action="store_true", default=False, help="if enabled, runs binary on global histogram")
parser.add_option("-q", "--queue", dest="queue", default="8nh", help="type of queue; if not specified, default is 8nh")
parser.add_option("-j","--job", dest="job", default="job", help="job name; if not specified, default is job")
parser.add_option("-n","--number", dest="number", default="30", help="number of files per job; if not specified, default is 30")
parser.add_option("-s","--storage", dest="storage", default="CASTOR", help="storage of files to process; if not specified, default is CASTOR")
parser.add_option("-p","--path", dest="path", default="/afs/cern.ch/na62/offline/lists/Data/Reco/2016/", help="path of files (on CASTOR/EOS) to be processed; if not specified, default is /afs/cern.ch/na62/offline/lists/Data/Reco/2016/")
parser.add_option("-o","--output", dest="output", default="./MyAnalyzerOutput", help="path of output files; if not specified, default is ./MyAnalyzerOutput")

(options, args) = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

arg.add_argument('args', nargs='*')

binary = os.path.abspath(sys.argv[1])
binary_name = os.path.split(sys.argv[1])[1]
list_of_runs = os.path.abspath(sys.argv[2])
launch = options.launch
add = options.add
run = options.run
queue = options.queue
job = options.job
number = options.number
storage = options.storage
path = options.path
output = "/".join(os.path.abspath(options.output).split('/')[:-1])
name = os.path.basename(os.path.normpath(options.output))

counter = 0;

# Function for job handling

def job_manager(binary, binary_name, path, run_name, number, job, queue):
    prefix = "xfarmer"
    list_of_files = os.listdir(path)
    files_in_dir = len(list_of_files)
    useless_files = 0
    for file in list_of_files:
        if run_name in file and "Electron" not in file and "3Tracks" not in file and storage.upper() in file:
            file_name = str(file)
        elif run_name not in file:
            useless_files += 1
    if files_in_dir == useless_files:
        print(str(run_name + " not found!"))
        shutil.rmtree(str(output + "/" + name + "/" + run_name))
        return
    os.system(str("grep -v ! " + path + file_name + "| split -a 3 -l " + number + " - " + prefix))
    current_dir = os.getcwd()
    farmer_in_run = os.listdir(current_dir)
    print(farmer_in_run)
    if len(farmer_in_run) == 0:
        print(str(run_name + " list has no items!"))
        shutil.rmtree(str(output + "/" + name + "/" + run_name))
        return
    if len(farmer_in_run) != 0:
        print(str(run_name + ": job submitted successfully!"))
    nj = 0
    os.chdir(current_dir)

    for item in farmer_in_run:
        if prefix not in item:
            continue
        elif prefix in item:
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
            with open("name.txt","w") as nn:
                nn.write(njfull)
            os.chdir(current_dir)
            nj += 1
    return

# Start

os.chdir(output)

if launch:
    print(str("[CheckJobs] About to execute binary " + binary_name + ". Output file is " + output + "/" + name))
if add:
    print("[CheckJobs] About to hadd all output histograms. Outputs will be in run directories")
if run:
    print(str("[CheckJobs] About to re-run " + binary_name + " over its own outputs. Outputs will be in run directories"))

# Parse list of runs

with open(list_of_runs) as f:
    for line in f:
        regex = re.match("(\d+)", line.rstrip())
        if regex:
            run_name = str("Run" + regex.group(1))
            if launch:                 # execute jobs
                if not os.path.exists(os.path.abspath(name)) and counter == 0:
                    os.makedirs(str(output + "/" + name))
                elif os.path.exists(os.path.abspath(name)) and counter == 0:
                    print("Directory for runs already exists! Remove it before executing this script")
                    sys.exit(1)
                os.chdir(str(output + "/" + name))
                os.makedirs(run_name)
                os.chdir(run_name)
                job_manager(binary, binary_name, path, run_name, number, job, queue)
                os.chdir(str(output + "/" + name))
            if add:                   # add histos in each run directory
                os.chdir(str(output + "/" + name))
                if os.path.exists(run_name):
                    os.chdir(run_name)
                if os.path.exists(str("Result.root")):
                    os.remove(str("Result.root"))
                if os.path.exists("outFile.root"):
                    os.remove("outFile.root")
                histos_to_add = " ".join(glob.glob('./*/outFile.root'))
                os.system(str("hadd Result.root " + histos_to_add))
                os.chdir(str(output + "/" + name))
        counter += 1

os.chdir(str(output + "/" + name))

if add:                              # add all histos
    print(str("[CheckJobs] About to hadd histograms from all run directories. Output will be in " + output + "/" + name))
    os.chdir(str(output + "/" + name))
    if os.path.exists(str("Result" + binary_name + ".root")):
        os.remove(str("Result" + binary_name + ".root"))
    histos_to_add = " ".join(glob.glob('./*/Result.root'))
    output_histo = str("Result" + binary_name + ".root")
    os.system(str("hadd " + output_histo + " " + histos_to_add))
if run:                              # re-run analyzer on final histo
    print(str("[CheckJobs] About to re-run " + binary_name + " over global output. Output will be in " + output + "/" + name))
    os.chdir(str(output + "/" + name))
    if not os.path.exists(str("Result" + binary_name + ".root")):
        print("Cannot run, input file does not exists!")
        sys.exit(1)
    if os.path.exists(str("Result" + binary_name + "2.root")):
        os.remove(str("Result" + binary_name + "2.root"))
    os.system(str(binary + " -i " + "Result" + binary_name + ".root -o Result" + binary_name + "2.root --histo"))
