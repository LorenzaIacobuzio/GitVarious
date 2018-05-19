#!/usr/bin/env bash

scp liacobuz@lxplus.cern.ch:~/NA62AnalysisTool/bla*.root liacobuz@lxplus.cern.ch:~/NA62AnalysisTool/minc*.root ~/Desktop
cd ~/Desktop/GitVarious/MyScripts/Others/AnalyzerPlots/
#root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/bla.root", "~/Desktop/minc.root")'
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/bla1.root", "~/Desktop/minc1.root")'
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/bla2.root", "~/Desktop/minc2.root")'
root -b -q 'Plots.C("~/Desktop/MCnote/images/Plots/", "~/Desktop/bla3.root", "~/Desktop/minc3.root")'
