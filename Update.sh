#!/usr/bin/env bash

scp -r liacobuz@lxplus.cern.ch:/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/Analyzers .
scp -r liacobuz@lxplus.cern.ch:/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/HaddScript MyScripts/Others/
scp liacobuz@lxplus.cern.ch:/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/CondorJobs/* MyScripts/Jobs/
git status
git add -A
git commit -m "update"
git push origin master
