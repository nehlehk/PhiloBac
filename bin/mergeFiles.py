#!/usr/bin/env python

import fileinput
import sys


inputFiles = []
for i in range(1,len(sys.argv)):
    inputFiles.append(sys.argv[i])


with open('AllOtherTrees.newick', 'w') as file:
    input_lines = fileinput.input(inputFiles)
    file.writelines(input_lines)


