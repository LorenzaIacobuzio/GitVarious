#!/bin/bash

cd /afs/cern.ch/work/l/liacobuz/private/condorOutput/HNL/outFiles

for i in `seq 1 10`;
	do
		hadd -f final$i.root $i/new*.root
	done

for i in `seq 1 10`;
        do
		hadd -f Final.root final*.root
