#!/bin/bash

cd /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool
source /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh
~/NA62AnalysisTool/MyAnalyzer -l ~/NA62AnalysisTool/List.txt -o /afs/cern.ch/work/l/liacobuz/private/HNLMCAnalysisAccYield.root
