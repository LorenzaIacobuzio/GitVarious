#!/bin/bash

## Functions

function findListFile() {

    local listFilePattern=""
    local swVersion=0
    local latestSWVersion=0
    local regex='[^0-9]'

    ## Clear some previous variable values
    listFile=""
    year=""
    isData=""
    version=""

    ## Check if we are working with data or MC
    if [[ ${run} =~ ${regex} ]]; then
	   isData=0
    else
	   isData=1
    fi

    echo "isData: ${isData}"
    if [ ${isData} -eq 1 ]; then

	   ## Check if we are running with filtered or unfiltered data
	   if [ "$filter" == "nofilter" ]; then
	       #listFilePattern="Run${run}Raw.EOSlist.r"
	       #listFilePattern="Run${run}Raw.EOSlist.r"
               #Modification for reprocessing unfiltered lists stored in my directory
	       listFilePattern="Run00${run}.EOSlist.${filter}.p"
	   else
	       #listFilePattern="Run${run}Raw.EOSlist.${filter}.r"
               #Modification for reprocessing > r1666
	       listFilePattern="Run00${run}.EOSlist.${filter}.p"
               #Modification for reprocessing stored in Karim's directory
	       #listFilePattern="Run00${run}.${filter}.r"

	   fi

        ## Determine from which year is the run
        if (( run >= 3015 && run <= 4174 )); then
            year=2015
        elif (( run >= 4560 && run <= 6912 )); then
            year=2016
        elif (( run >= 7281 && run <= 8389 )); then
            year=2017
        elif (( run >= 8390 )); then
            year=2018
        else
            year=0
        fi

	if [[ ${sample} == *"2016A"* ]]; then
	    version="v1.0.5"
	elif [[ ${sample} == *"2017A"* ]]; then
	    version="v1.0.2"
	elif [[ ${sample} == *"2017B"* ]]; then
	    version="v1.0.3"
	elif [[ ${sample} == *"2017C"* ]]; then
	    version="v1.0.3"
	elif [[ ${sample} == *"2017D"* ]]; then
	    version="v1.0.3"
	elif [[ ${sample} == *"2018A"* ]]; then
	    version="v1.1.2"
	elif [[ ${sample} == *"2018B"* ]]; then
	    version="v1.1.3"
	elif [[ ${sample} == *"2018E"* ]]; then
	    version="v1.0.4"
	elif [[ ${sample} == *"2018H"* ]]; then
	    version="v1.0.5"
	fi

	echo "sample and version: " ${sample} ${version}

        echo "year: ${year}"
        if [ ${year} -eq 0 ]; then
            echo "Run: ${run} could not be assigned a year..."
            return
        else
	    listsLocation=${dataListsLocation}/${sample}-${version}
        fi

        echo "list location: ${listsLocation}"

        local listCandidates=( $(ls ${listsLocation} | grep ${listFilePattern}) )
        for i in "${listCandidates[@]}"
        do
            if [ "$filter" == "nofilter" ]; then
                #swVersion=${i#*.r}
		# Modification for reprocessing unfiltered stored in my directory
                swVersion=${i#*.p}
		echo "SW version: ${swVersion}"
            else
                # TODO: Add check for filter revision
                #swVersion=${i#*.r}
                #swVersion=${swVersion%%.r*}
                # Modification for reprocessing > r1666
                swVersion=${i#*.p}
                # Modification for reprocessing stored in Karim's directory
                #swVersion=${i#*.r}
		echo "SW version: ${swVersion}"
            fi

            # Find the latest available sw revision
            #if (( swVersion > latestSWVersion )); then
            latestSWVersion=$swVersion
            listFile=${listsLocation}/${i}
	    echo $listFile
            #fi
        done
    else
	   listFile=${mcListsLocation}/${run}.list
    fi
}

## Check if an input was supplied
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters. Exiting..."
    exit 0
fi
configFile="$1"

## Check if config file exists
if [ ! -f ${configFile} ]; then
    echo "Config file not found. Exiting..."
    exit 0
fi

echo "Parameters from config file:"

## Parse config file
while read -r configLine || [[ -n "$configLine" ]]; do

    # Skip lines begining with #
    if [[ ${configLine:0:1} == "#" ]]; then
	continue
    fi

    # Find data lists location
    dataListsLocationPrefix="dataListsLocation= "
    if [[ ${configLine} == ${dataListsLocationPrefix}* ]]; then
	dataListsLocation=${configLine#$dataListsLocationPrefix}
	echo "Data lists dir:" ${dataListsLocation}
    fi

    # Find MC lists location
    mcListsLocationPrefix="mcListsLocation= "
    if [[ ${configLine} == ${mcListsLocationPrefix}* ]]; then
	mcListsLocation=${configLine#$mcListsLocationPrefix}
	echo "MC lists dir:" ${mcListsLocation}
    fi

    # Find analysis script name
    analysisScriptPrefix="analysisScript= "
    if [[ ${configLine} == ${analysisScriptPrefix}* ]]; then
	analysisScript=${configLine#$analysisScriptPrefix}
	echo "Script to run:" ${analysisScript}
    fi

    # Find executable name
    executableNamePrefix="executableName= "
    if [[ ${configLine} == ${executableNamePrefix}* ]]; then
	executableName=${configLine#$executableNamePrefix}
	echo "Name of executable:" ${executableName}
    fi

    # Find executable file
    executableFilePrefix="executableFile= "
    if [[ ${configLine} == ${executableFilePrefix}* ]]; then
	executableFile=${configLine#$executableFilePrefix}
	echo "Executable:" ${executableFile}
    fi

    # Find runs
    runsPrefix="runs= "
    if [[ ${configLine} == ${runsPrefix}* ]]; then
	runs=${configLine#$runsPrefix}
	echo "Runs:" ${runs}
    fi

    # Find filter
    filterPrefix="filter= "
    if [[ ${configLine} == ${filterPrefix}* ]]; then
	filter=${configLine#$filterPrefix}
	echo "Filter:" ${filter}
    fi

    # Find file multiplier
    filesPerJobPrefix="filesPerJob= "
    if [[ ${configLine} == ${filesPerJobPrefix}* ]]; then
	filesPerJob=${configLine#$filesPerJobPrefix}
	echo "Files per job:" ${filesPerJob}
    fi

    # Find job flavour
    jobFlavourPrefix="jobFlavour= "
    if [[ ${configLine} == ${jobFlavourPrefix}* ]]; then
    jobFlavour=${configLine#$jobFlavourPrefix}
    echo "Job flavour:" ${jobFlavour}
    fi

    # Find number of CPUs
    numberOfCPUsPrefix="numberOfCPUs= "
    if [[ ${configLine} == ${numberOfCPUsPrefix}* ]]; then
    numberOfCPUs=${configLine#$numberOfCPUsPrefix}
    echo "Number of CPUs:" ${numberOfCPUs}
    fi

    # Find output directory
    outputDirPrefix="outputDir= "
    if [[ ${configLine} == ${outputDirPrefix}* ]]; then
	outputDir=${configLine#$outputDirPrefix}
	echo "Output dir:" ${outputDir}
    fi

    # Find sample
    samplePrefix="outputDir= "
    if [[ ${configLine} == ${samplePrefix}* ]]; then
	sample=${configLine: -6:5}
	echo "Sample:" ${sample}
    fi

    # Find analysis directory
    analysisDirPrefix="analysisDir= "
    if [[ ${configLine} == ${analysisDirPrefix}* ]]; then
	analysisDir=${configLine#$analysisDirPrefix}
	echo "Analysis dir:" ${analysisDir}
    fi

    # Find my lists directory
    myListDirPrefix="myListDir= "
    if [[ ${configLine} == ${myListDirPrefix}* ]]; then
	myListDir=${configLine#$myListDirPrefix}
	echo "My list dir:" ${myListDir}
    fi

    # Find executable options
#    myOptionsPrefix="myOptions= "
#    if [[ ${configLine} == ${myOptionsPrefix}* ]]; then
#	myOptions=${configLine#$myOptionsPrefix}
#	echo "My options:" ${myOptions}
#    fi

   # Find executable ignore
    myIgnorePrefix="myIgnore= "
    if [[ ${configLine} == ${myIgnorePrefix}* ]]; then
        myIgnore=${configLine#$myIgnorePrefix}
        echo "My ignore:" ${myIgnore}
    fi

done < <(cat ${configFile})

## Keep a running counter of total jobs submitted
totalJobs=0

## Parse list of runs, and submit jobs
for run in ${runs}; do
    echo ""
    echo "Processing: ${run}"

    ## Check if we are working with a run/MC, or a list of files
    if [ "${run: -5}" == ".list" ]; then
        echo "Working with a list of files, not a run"
        filter=nofilter
        isData=1
        listFile=${myListDir}/${run}
        run="${run%.*}"

        ## Check if the list of files physically exists
        if [ ! -f ${listFile} ]; then
            echo "List of files, ${listFile}, not found. Skipping..."
            continue
        fi
    else
        ## Function to find the relevant list file, based on run number and filter, or MC sample.
        findListFile "$sample"

        ## Check if a list was found
        if [ -z "${listFile}" ]; then
            echo "Could not find list file for run: ${run}, and filter: ${filter}. Skipping ..."
            continue
        fi
        echo "Year of run: ${year}"
        echo "List file: ${listFile}"
    fi
    # Append filter to dir name, if it's data and there is a filter
    if [ ${isData} -eq 1 ] && [ "${filter}" != "nofilter" ]; then
        runOuputDir=${outputDir}/${run}_${filter}
    else
        runOuputDir=${outputDir}/${run}
    fi
    echo "Run output dir: ${runOuputDir}"

    # Check that the run output dir doesn't exist already
    while [ -d ${runOuputDir} ]; do
        read -p "The above run output dir already exists. Would you like to overwrite it? (y/n)" yn
        case $yn in

            [Yy]* ) echo "OK! I will overwrite it."
                    rm -r ${runOuputDir}
		    rm /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/NEventsFile.txt
		    rm /afs/cern.ch/user/l/liacobuz/NA62AnalysisTool/SumGoodFile.txt
                    ;;
            [Nn]* ) read -p "OK. Please enter a new run dir: " newRunOutputDir
                    runOuputDir=${outputDir}/${newRunOutputDir}
                    echo "New run output dir: ${runOuputDir}"
                    ;;
            * )     echo "Please answer y or n."
                    ;;
        esac
    done

    # Create the directories for this analysis
    cutListsDir=${runOuputDir}/cutLists
    mkdir -p ${cutListsDir}

    # Create directories for condor files
    mkdir -p ${runOuputDir}/error
    mkdir -p ${runOuputDir}/log
    mkdir -p ${runOuputDir}/output

    inputDir=${runOuputDir}/input
    mkdir -p ${inputDir}

    # Full run check, at the moment forcing full run...
    firstLine=0
    if [ ${firstLine} -eq 0 ]; then
	   firstLine=1
	   numberOfFiles=$(wc -l < ${listFile})
	   echo "Number of files in the list: ${numberOfFiles}"
	   lastLine=$(((${numberOfFiles}+(${filesPerJob}-1))/${filesPerJob}))
	   echo "I will use the full list"
    else
	   echo "First counter: " ${firstLine}
	   echo "Last counter: " ${lastLine}
    fi
    echo "There will be ${lastLine} jobs for this run"

    ## Create individual list files, and analysis script input file
    analysisInputFile=analysisInput_${firstLine}_${lastLine}.txt

    for i in `seq ${firstLine} ${lastLine}`;
    do
	   jobNumber=$((${i}*${filesPerJob}))
	   jobListFile=${cutListsDir}/${run}_${jobNumber}.list

        if [ "${jobNumber}" -lt "${numberOfFiles}" ];  then
	       head -${jobNumber} ${listFile} | tail -${filesPerJob} > ${jobListFile}
        else
            finalList=$((${filesPerJob} - (${jobNumber}-${numberOfFiles}) ))
            head -${jobNumber} ${listFile} | tail -${finalList} > ${jobListFile}
        fi

	   #analysisJobInput="${executableFile} ${isData} ${run} ${jobNumber} ${analysisDir} ${jobListFile} ${myOptions}"
analysisJobInput="${executableFile} ${isData} ${run} ${jobNumber} ${analysisDir} ${jobListFile} ${myIgnore}"
	   echo ${analysisJobInput} >> ${inputDir}/${analysisInputFile}
    done

    ## Create condor job submission file
    condorSubmissionFile=${runOuputDir}/${executableName}_${firstLine}_${lastLine}.sub
    executableLine="executable              = ${analysisScript}"
    outputLine="output                  = output/\$(ClusterId).\$(ProcId).out"
    errorLine="error                   = error/\$(ClusterId).\$(ProcId).err"
    logLine="log                     = log/\$(ClusterId).log"
    initialDirLine="initialdir              = ${runOuputDir}"
    outputFileLine="transfer_output_files   = outFiles, usedLists"
    cpuLine="RequestCpus             = ${numberOfCPUs}"
    jobFlavourLine="+JobFlavour             = \"${jobFlavour}\""
    queueLine="queue arguments from ${inputDir}/${analysisInputFile}"

    echo "${executableLine}" >> ${condorSubmissionFile}
    echo "${outputLine}" >> ${condorSubmissionFile}
    echo "${errorLine}" >> ${condorSubmissionFile}
    echo "${logLine}" >> ${condorSubmissionFile}
    echo "${initialDirLine}" >> ${condorSubmissionFile}
    echo "${outputFileLine}" >> ${condorSubmissionFile}
    echo "${cpuLine}" >> ${condorSubmissionFile}
    echo "${jobFlavourLine}" >> ${condorSubmissionFile}
    echo "" >> ${condorSubmissionFile}
    echo "${queueLine}" >> ${condorSubmissionFile}

    ## Submit the job
    condor_submit ${condorSubmissionFile}

    ((totalJobs = totalJobs + lastLine))
done

echo ""
echo "${totalJobs} jobs submitted in total"
