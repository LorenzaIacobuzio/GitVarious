#!/usr/bin/env bash

scp liacobuz@lxplus.cern.ch:~/NA62AnalysisTool/*.root ~/Desktop
cd ~/Desktop/GitVarious/MyScripts/Others/AnalyzerPlots/
root -b -q 'Plots.C("", "~/Desktop/1.root", false)'
root -b -q 'Plots.C("", "~/Desktop/2.root", false)'
root -b -q 'Plots.C("", "~/Desktop/3.root", false)'
root -b -q 'Plots.C("", "~/Desktop/ToyMC.root", true)'
