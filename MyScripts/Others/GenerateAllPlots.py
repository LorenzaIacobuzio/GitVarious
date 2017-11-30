#!/usr/bin/env python

import os
import sys

os.system("/li/home/GitVarious/MyScripts/Others/root -l")

for i in range(0, 4):
    for j in range(1, 4):
        os.system(".L GeneralPlots.C")
        os.system("GeneralPlots({0}, {1}".format(i, j))
