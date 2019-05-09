#!/usr/bin/env bash

scp -r $LX/Analyzers .
scp -r $LX/HaddScript MyScripts/Others/
scp $LX/CondorJobs/* MyScripts/Jobs/
git status
git add -A
git commit -m "update"
git push origin master
