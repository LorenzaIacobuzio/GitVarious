#!/bin/env python

import os
import sys
import subprocess

def createJobs(njobs):
	for i in xrange(1,njobs+1):
		subprocess.call(["bsub", "-q", "8nm", "-oo", "$WORK/log{0}".format(i),"""#!/bin/bash   cd /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/
source /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/scripts/env.sh
cd /afs/cern.ch/user/l/liacobuz/na62fw/NA62Reconstruction/
source /afs/cern.ch/user/l/liacobuz/na62fw/NA62Reconstruction/scripts/env.sh
cd /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool
source /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh
/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/MyAnalyzer -i $WORK/SingleHistos/Reco/HNLMCReco_{0}.root -o $WORK/HNLMCAnalysis_{0}.root""".format(i)])
		

if __name__=="__main__":
	njobs = int(sys.argv[1])
	createJobs(njobs)


