#!/bin/bash

cd /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool
source /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh
~/NA62AnalysisTool/MyAnalyzer -i $WORK/ImportantHistos/Ordinary/HNLMCReco.root -o $WORK/HNLMCAnalysis.root
