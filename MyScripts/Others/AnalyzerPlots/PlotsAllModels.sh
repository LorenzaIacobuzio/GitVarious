#!/usr/bin/env bash

#scp liacobuz@lxplus.cern.ch:~/NA62AnalysisTool/First*.root liacobuz@lxplus.cern.ch:~/NA62AnalysisTool/Second*.root ~/Desktop
cd ~/Desktop/GitVarious/MyScripts/Others/AnalyzerPlots/
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/First1.root", "~/Desktop/Second1.root")'
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/First2.root", "~/Desktop/Second2.root")'
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/First3.root", "~/Desktop/Second3.root")'
