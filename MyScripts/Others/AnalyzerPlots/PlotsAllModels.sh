#!/usr/bin/env bash

#scp liacobuz@lxplus.cern.ch:/afs/cern.ch/work/l/liacobuz/private/1.root liacobuz@lxplus.cern.ch:/afs/cern.ch/work/l/liacobuz/private/2.root liacobuz@lxplus.cern.ch:/afs/cern.ch/work/l/liacobuz/private/3.root liacobuz@lxplus.cern.ch:/afs/cern.ch/work/l/liacobuz/private/Sublist-POT.root liacobuz@lxplus.cern.ch:/afs/cern.ch/work/l/liacobuz/private/Toy.root ~/Desktop/NewHistos
cd ~/GitVarious/MyScripts/Others/AnalyzerPlots/
#root -b -q 'Plots.C("", "~/Desktop/NewHistos/1.root", false, false)'
root -b -q 'Plots.C("", "~/Desktop/NewHistos/2.root", false, false)'
#root -b -q 'Plots.C("", "~/Desktop/NewHistos/3.root", false, false)'
#root -b -q 'Plots.C("", "~/Desktop/NewHistos/POT.root", false, false)'
#root -b -q 'Plots.C("", "~/Desktop/NewHistos/Toy.root", true, false)'
#root -b -q 'Plots.C("", "~/Desktop/NewHistos/Data.root", false, true)'
