#!/usr/bin/env python

#########################################
# Author: Jianhong Ou
# This source code is licensed under the MIT license
#########################################

from __future__ import print_function
import argparse
import os
import sys
import io
import errno
import glob


def wrap_ini(matrix, tads):
    contents = """[x-axis]
where = top

[hic matrix]
file = {0}
title = Hi-C data
# depth is the maximum distance plotted in bp. In Hi-C tracks
# the height of the track is calculated based on the depth such
# that the matrix does not look deformed
depth = 300000
transform = log1p
file_type = hic_matrix
height = 9

[tads]
file = {1}
file_type = domains
border_color = black
color = none
# the tads are overlay over the hic-matrix
# the share-y options sets the y-axis to be shared
# between the Hi-C matrix and the TADs.
overlay_previous = share-y

    """.format(
        matrix, tads
    )
    return contents


def parse_args(args=None):
    Description = "Create HiCAR TADS track ini file"
    Epilog = """Example usage: python create_hicar_track_ini.py <COOLER_OUT> <TAD_BED>"""
    argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    ## REQUIRED PARAMETERS
    argParser.add_argument("--matrix", help="h5 matrix")
    argParser.add_argument("--tads", help="tad bed file")
    argParser.add_argument("--out", help="output ini filename")
    return argParser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    ini = wrap_ini(args.matrix, args.tads)
    out = open(args.out, "a")
    out.write(ini)
    out.close()


if __name__ == "__main__":
    sys.exit(main())
