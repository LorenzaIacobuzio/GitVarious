#!/usr/bin/env bash

scp -r $LX/Analyzers .
scp $LX/RunJobs.py MyScripts/Jobs/Various
git status
git add -A
git commit -m "update"
git push origin master
