#!/usr/bin/env bash

cd ~/NA62AnalysisTool
source ~/NA62AnalysisTool/scripts/env.sh
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/Reco.root -o $11.root -p "HeavyNeutrinoScan:UeSquaredRatio=52;UmuSquaredRatio=1;UtauSquaredRatio=1&HeavyNeutrino:UeSquaredRatio=52;UmuSquaredRatio=1;UtauSquaredRatio=1"
./bin-slc6/MyAnalyzer -i $11.root -o $21.root --histo
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/Reco.root -o $12.root -p "HeavyNeutrinoScan:UeSquaredRatio=1;UmuSquaredRatio=16;UtauSquaredRatio=3.8&HeavyNeutrino:UeSquaredRatio=1;UmuSquaredRatio=16;UtauSquaredRatio=3.8"
./bin-slc6/MyAnalyzer -i $12.root -o $22.root --histo
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/Reco.root -o $13.root -p "HeavyNeutrinoScan:UeSquaredRatio=0.061;UmuSquaredRatio=1;UtauSquaredRatio=4.3&HeavyNeutrino:UeSquaredRatio=0.061;UmuSquaredRatio=1;UtauSquaredRatio=4.3"
./bin-slc6/MyAnalyzer -i $13.root -o $23.root --histo
