#!/bin/bash

cd /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool
source /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh
~/NA62AnalysisTool/MyAnalyzer -i $WORK/ImportantHistos/.root -o $WORK/ImportantHistos/.root
