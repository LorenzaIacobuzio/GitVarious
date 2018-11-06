#!/bin/bash



run=$1
type=$2
roottools=/afs/cern.ch/work/p/piandani/NA62/na62fw/Analysis/jobs/mergedfiles/copyFiles
analyzer=$3

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