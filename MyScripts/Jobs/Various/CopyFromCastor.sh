#!/bin/bash
source /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/scripts/env.sh
#files=`cat /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/HeavyNeutrino/HNList/k3pi.list`
files=`cat /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/eos_file.txt`
for f in $files; do
#	$XROOTD/bin/xrdcp xroot://castorpublic.cern.ch/$f .
	$XROOTD/bin/xrdcp xroot://eosna62.cern.ch//eos/na62/user/l/liacobuz
done

ls *.root > localList.lst
/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/MyAnalyzer -l localList.lst -o out.root 
#cp out.root /afs/cern.ch/work/l/liacobuz/private/HNAnOnK3pi.root
