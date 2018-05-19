#!/usr/bin/env bash

cd ~/NA62AnalysisTool
source ~/NA62AnalysisTool/scripts/env.sh
NA62AnalysisBuilder.py configfiles/MyConfig -j16
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/2.root -o bla1.root -p "HeavyNeutrinoScan:UeSquaredRatio=52;UmuSquaredRatio=1;UtauSquaredRatio=1&HeavyNeutrino:UeSquaredRatio=52;UmuSquaredRatio=1;UtauSquaredRatio=1" -n 100
./bin-slc6/MyAnalyzer -i bla1.root -o minc1.root --histo
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/2.root -o bla2.root -p "HeavyNeutrinoScan:UeSquaredRatio=1;UmuSquaredRatio=16;UtauSquaredRatio=3.8&HeavyNeutrino:UeSquaredRatio=1;UmuSquaredRatio=16;UtauSquaredRatio=3.8" -n 100
./bin-slc6/MyAnalyzer -i bla2.root -o minc2.root --histo
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/2.root -o bla3.root -p "HeavyNeutrinoScan:UeSquaredRatio=0.061;UmuSquaredRatio=1;UtauSquaredRatio=4.3&HeavyNeutrino:UeSquaredRatio=0.061;UmuSquaredRatio=1;UtauSquaredRatio=4.3" -n 100
./bin-slc6/MyAnalyzer -i bla3.root -o minc3.root --histo
