#!/usr/bin/env bash

scp liacobuz@lxplus.cern.ch:~/NA62AnalysisTool/$1*.root liacobuz@lxplus.cern.ch:~/NA62AnalysisTool/$2*.root ~/Desktop
cd ~/Desktop/GitVarious/MyScripts/Others/AnalyzerPlots/
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/$11.root", "~/Desktop/$21.root")'
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/$12.root", "~/Desktop/$22.root")'
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/$13.root", "~/Desktop/$23.root")'
