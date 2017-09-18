#!/usr/bin/env python

import re
import datetime
import collections

def CreateCSV(detectorName):
    boards = dict()
    with open("GlobalBoardReload.txt", "r") as f:
        for line in f:
            regex = "([\d.]+)\s([\d+:.]+)[\s\w]+TELL_" + detectorName + "(_(\d+))?"
            reg = re.search(regex, line, re.IGNORECASE)
            if reg:
                if reg.group(4) is not None:
                    boardnumber = int(reg.group(4))
                else:
                    boardnumber = 0
                y,m,d = [int(x) for x in reg.group(1).split('.')]
                if y == 2017 and m >= 5:
                    weeknumber = datetime.date(y,m,d).isocalendar()[1] - 19
                    if weeknumber not in boards:
                        boards[weeknumber] = {}
                        if boardnumber not in boards[weeknumber]:
                            boards[weeknumber][boardnumber] = 1
                        else:
                            boards[weeknumber][boardnumber] += 1
                    else:
                        if boardnumber not in boards[weeknumber]:
                            boards[weeknumber][boardnumber] = 1
                        else:
                            boards[weeknumber][boardnumber] += 1

    orderedboards = collections.OrderedDict(sorted(boards.items()))

    with open("BoardReloads" + str(detectorName) + ".csv", "w") as csv:
        for key in orderedboards.iterkeys():
            csv.write(str(key))
            for i in range(6):
                if type(orderedboards[key].get(i)) is not int:
                    orderedboards[key][i] = 0
                csv.write("," + str(orderedboards[key].get(i)))
            csv.write("\n")

    return

if __name__ == "__main__":
    detectors = ["KTAG", "CHANTI", "LAV1",  "LAV2", "LAV3", "LAV4", "LAV5", "LAV6", "LAV7", "LAV8", "LAV9", "LAV10", "LAV11", "LAV12", "IRC_SAC", "MUV3", "HASC", "RICH"]
    for item in detectors:
        CreateCSV(item)
