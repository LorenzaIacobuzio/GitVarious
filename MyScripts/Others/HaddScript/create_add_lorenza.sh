#!/bin/bash

# to run: ./exe run type analyzer

run=$1
type=$2 # data/mc (doesn't really matter)
roottools=/afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/HaddScript/copyFiles
analyzer=$3 # analyzer name

incuttedfile=root$run.list
if  [[ -e $incuttedfile ]];
then
    echo "input file exist... removing it"
    rm -f $incuttedfile
fi

ls $run/ | grep *.root > root$run.list
ofile=Addroot$run

echo "Adding file for run "$run 

if  [[ -e $ofile ]];
then
    echo "outfile exist... removing it"
    rm -f $ofile
fi

echo -n "hadd test.root ">$ofile
counter=0
while read line ; do
    echo "processing file "$counter
    counter=$(($counter+1))
    $roottools $run/$line $run/$type\_$line $analyzer
    echo -n $run/$type\_$line" " >>$ofile

done < "$incuttedfile"

chmod 777 $ofile
./$ofile