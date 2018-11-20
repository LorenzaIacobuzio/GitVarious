#!/usr/bin/env bash

scp -r $LX/Analyzers .
git status
git add -A
git commit -m "update"
git push origin master
