#!/bin/bash

cd /afs/cern.ch/work/l/liacobuz/private/$1/HNL/outFiles

i=1;while read l;do mkdir $i;mv $l $((i++));done< <(ls|xargs -n100)

