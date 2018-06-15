#!/usr/bin/env bash

#scp liacobuz@lxplus.cern.ch:/afs/cern.ch/work/l/liacobuz/private/*.root ~/Desktop
cd ~/Desktop/GitVarious/MyScripts/Others/AnalyzerPlots/
root -b -q 'Plots.C("", "~/Desktop/Working/1.root", false)'
root -b -q 'Plots.C("", "~/Desktop/Working/2.root", false)'
root -b -q 'Plots.C("", "~/Desktop/Working/3.root", false)'
root -b -q 'Plots.C("", "~/Desktop/Working/ToyMC.root", true)'
