#!/bin/bash

## Input check
if [ $# -lt 6 ]; then
   echo "$# is a illegal number of parameters!"
   exit 0
fi

executableFile="$1"
echo "Executable file: ${executableFile}"

isData="$2"
echo "Working on data? ${isData}"

run="$3"
echo "Run number ${run}"

jobNumber="$4"
echo "Job number: ${jobNumber}"

analysisDir="$5"
echo "NA62 analysis directory: ${analysisDir}"

jobListFile="$6"
echo "List file: ${jobListFile}"

myOptions="$7 $8"
echo "Executable options: ${myOptions}"

## CASTOR
castorPrefix="root://castorpublic.cern.ch/"

## Local directories and files
mkdir inFiles
mkdir outFiles
mkdir usedLists
usedList="usedLists/${run}_${jobNumber}.list"
outputFile="outputFile.txt"

#Copy files over from eos or castor
while read -r line || [[ -n "$line" ]]; do
	fileCopiedOK=false
	copyAttempts=0
	fileName=${line##*/}
	if [[ ${line} == *"castor"* ]]; then
		echo "Working with files on CASTOR."
		while [ ${fileCopiedOK} == false ]; do
			## Don't try to copy the file indefinitely
			if [ ${copyAttempts} -ge 5 ]; then
				echo "Giving up copying ${fileName} after ${copyAttempts} attempts"
				break
			fi

			## Attempt to copy the file...
			xrdcp ${castorPrefix}${line} inFiles/${fileName} > ${outputFile} 2>&1
			cat ${outputFile}
			((copyAttempts++))
			if [ -f "inFiles/${fileName}" ]; then
				echo "${fileName} exists"
				## Check size
				castorFileSize=$(nsls -l ${line} | awk '{print $5}')
				localFileSize="$(wc -c <"inFiles/${fileName}")"
				if [ ${castorFileSize} -eq ${localFileSize} ]; then
					echo "${fileName} sizes match. File size: ${localFileSize}"
					echo "${fileName} copied OK"
					echo "inFiles/${line##*/}" >> ${usedList}
					fileCopiedOK=true
				else
					echo "${fileName} sizes do not match."
					echo "Source file size: ${castorFileSize}"
					echo "Copied file size: ${localFileSize}"

				fi
			else
				echo "${fileName} does not exist."
			fi
		done
	else
		echo "Working with files on EOS."
		while [ ${fileCopiedOK} == false ]; do
			## Don't try to copy the file indefinitely
			if [ ${copyAttempts} -ge 5 ]; then
				echo "Giving up copying ${fileName} after ${copyAttempts} attempts"
				break
			fi

			## Attempt to copy the file...
			cp ${line} inFiles/${fileName} > ${outputFile} 2>&1
			cat ${outputFile}
			((copyAttempts++))
			if [ -f "inFiles/${fileName}" ]; then
				echo "${fileName} exists"
				#eosFileSize=$(fileinfo ${line}| grep "Size: ")
				#eosFileSize=${eosFileSize#*'Size: '}
				localFileSize="$(wc -c <"inFiles/${fileName}")"
				#if [ ${eosFileSize} -eq ${localFileSize} ]; then
				#	echo "${fileName} sizes match. File size: ${localFileSize}"
				#	echo "${fileName} copied OK"
				echo "inFiles/${line##*/}" >> ${usedList}
				fileCopiedOK=true
				#else
				#	echo "${fileName} sizes do not match."
				#	echo "Source file size: ${eosFileSize}"
				#	echo "Copied file size: ${localFileSize}"
				#fi
			else
				echo "${fileName} does not exist."
			fi
		done
	fi
done < <(cat ${jobListFile})

numberOFiles=$(ls inFiles | wc -l)

if (( $numberOFiles == 0)); then
	echo "No files to process. Exiting..."
	exit 0
fi

## Run the analysis
source ${analysisDir}/scripts/env.sh

outFile=outFiles/${run}_${jobNumber}.root

#${executableFile} -l ${usedList} -p "ThreeTrackAnalysis:Reference=3" -o ${outFile} > ${outputFile} 2>&1
#${executableFile} -l ${usedList} -o ${outFile} > ${outputFile} 2>&1
${executableFile} -l ${usedList} -e 2 -o ${outFile} ${myOptions} > ${outputFile} 2>&1
cat ${outputFile}
cat NEventsFile.txt >> /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/NEventsFile.txt
cat SumGoodFile.txt >> /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/SumGoodFile.txt