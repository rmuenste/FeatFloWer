#/usr/bin/env python
# vim: set filetype=python
"""
A module for extracting the cmake build ids 
"""

import re
import os

#===============================================================================
#                        Function main
#===============================================================================
def main():

    with open("../Feat_FloWer/cmake/modules/SetFlagsForID.cmake") as f:

        for line in f:

            if re.search(r"Q2P1_BUILD_ID STREQUAL", line):
                words = line.strip().split(" ")
                ids = words[2].split('"')
                print("Build ID: " + ids[1])



#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == '__main__':
    main()