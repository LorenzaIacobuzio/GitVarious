#!/bin/bash

cd $1
ndir=`ls -lR | grep \d | wc -l`

for i in `seq 1 $ndir`;
	do
		hadd -f final$i.root $i/new*.root
	done

hadd -f Final.root final*.root
