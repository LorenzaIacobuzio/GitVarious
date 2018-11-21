#!/usr/bin/env bash

scp -r $LX/Analyzers .
scp $LX/Jobs/* MyScripts/Jobs
scp $LX/HaddScripts/* MyScripts/Others/HaddScript
git status
git add -A
git commit -m "update"
git push origin master
