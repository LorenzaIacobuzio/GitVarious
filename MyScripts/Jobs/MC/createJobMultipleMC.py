#!/bin/env python

import subprocess
import sys
from Cheetah.Template import Template

def createJobs(njobs, macro):
	tp = Template(file=macro)
	for i in xrange(1,njobs+1):
		tp.seed = i
		tp.n = i
	
		with open("jobMacros/mac{0}.mac".format(i), "w") as fd:
			fd.write(str(tp))
		subprocess.call(["bsub", "-q", "8nh", "-oo", "/afs/cern.ch/work/l/liacobuz/private/log{0}".format(i),"""#!/bin/bash

cd /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/
source /afs/cern.ch/user/l/liacobuz/na62fw/NA62MC/scripts/env.sh
NA62MC jobMacros/mac{0}.mac""".format(i)])
		
	

if __name__=="__main__":
	njobs = int(sys.argv[1])
	mac = sys.argv[2]
	createJobs(njobs, mac)
