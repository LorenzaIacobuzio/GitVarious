#!/bin/env python

import subprocess
import sys

def createJobs(njobs, nevt):
	for i in xrange(0,njobs):
		start = i*int(nevt)
		subprocess.call(["bsub", "-q", "1nh", "-oo", "log{0}".format(i),"""#!/bin/bash

cd /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/
source /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/scripts/env.sh
cd /afs/cern.ch/user/l/liacobuz/na62fw/NA62Reconstruction/
source /afs/cern.ch/user/l/liacobuz/na62fw/NA62Reconstruction/scripts/env.sh
./bin/NA62Reco -i $WORK/ImportantHistos/HNLMC.root -o $WORK/ImportantHistos/HNLMCReco_{0}.root -c config/NA62Reconstruction.2017.conf -j {1} -n {2}""".format(i, start, nevt)])
		
	

if __name__=="__main__":
	njobs = int(sys.argv[1])
	nevt = sys.argv[2]
	createJobs(njobs, nevt)
