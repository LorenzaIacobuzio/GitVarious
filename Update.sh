#!/usr/bin/env bash

scp -r $LX/Analyzers .
scp -r $LX/HaddScript MyScripts/Others/HaddScript
scp -r $LX/Jobs MyScripts/Jobs/Condor
git status
git add -A
git commit -m "update"
git push origin master
