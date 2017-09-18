#!/bin/env python

import subprocess
import sys

def createJobs(njobs):
	for i in xrange(0,njobs):
		subprocess.call(["bsub", "-q", "1nh", "-oo", "/afs/cern.ch/work/l/liacobuz/private/log{0}".format(i),"""#!/bin/bash

cd /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/
source /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/scripts/env.sh
cd /afs/cern.ch/user/l/liacobuz/na62fw/NA62Reconstruction/
source /afs/cern.ch/user/l/liacobuz/na62fw/NA62Reconstruction/scripts/env.sh
./bin/NA62Reco -i $WORK/SingleHistos/MC/HNLMC_{0}.root -o $WORK/SingleHistos/Reco/HNLMCReco_{0}.root -c config/NA62Reconstruction.2017.conf""".format(i)])
		
	

if __name__=="__main__":
	njobs = int(sys.argv[1])
	createJobs(njobs)
