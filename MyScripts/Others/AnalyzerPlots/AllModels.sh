#!/usr/bin/env bash

cd ~/NA62AnalysisTool
source ~/NA62AnalysisTool/scripts/env.sh
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/Reco.root -o /afs/cern.ch/work/l/liacobuz/private/1.root -p "HeavyNeutrinoScan:UeSquaredRatio=52;UmuSquaredRatio=1;UtauSquaredRatio=1&HeavyNeutrino:UeSquaredRatio=52;UmuSquaredRatio=1;UtauSquaredRatio=1"
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/Reco.root -o /afs/cern.ch/work/l/liacobuz/private/2.root -p "HeavyNeutrinoScan:UeSquaredRatio=1;UmuSquaredRatio=16;UtauSquaredRatio=3.8&HeavyNeutrino:UeSquaredRatio=1;UmuSquaredRatio=16;UtauSquaredRatio=3.8"
./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/Reco.root -o /afs/cern.ch/work/l/liacobuz/private/3.root -p "HeavyNeutrinoScan:UeSquaredRatio=0.061;UmuSquaredRatio=1;UtauSquaredRatio=4.3&HeavyNeutrino:UeSquaredRatio=0.061;UmuSquaredRatio=1;UtauSquaredRatio=4.3"
#./bin-slc6/MyAnalyzer -i /afs/cern.ch/work/l/liacobuz/private/ImportantHistos/Reco/ToyMC.root -o /afs/cern.ch/work/l/liacobuz/private/ToyMC.root
