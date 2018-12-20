#!/bin/bash

source /afs/cern.ch/user/l/liacobuz/na62fw/NA62Reconstruction/scripts/env.sh
source /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh
cp /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/HaddScript/Hadd.py .
echo $1 $2 $3
./Hadd.py $1 $2 $3
