MC:

Use NA62AnalysisTool/CondorJobs to run HeavyNeutrino+HeavyNeutrinoScan analyzer on reconstructed MC files (set all parameters)
Use NA62AnalysisTool/HaddScript/MCHadd/MCHadd.py to copyFiles and keep only HeavyNeutrinoScan directory in root files
Run analyzer in --histo mode on hadded final root file

Data:

Use NA62AnalysisTool/CondorJobs to run HeavyNeutrino+HeavyNeutrinoScan analyzer on data
Use NA62AnalysisTool/HaddScript/DataHadd/DataHadd.py to have a final hadded root file
Run analyzer in --histo mode on hadded final root file
