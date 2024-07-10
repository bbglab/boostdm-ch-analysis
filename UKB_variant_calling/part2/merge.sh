#!/bin/sh

#first, transfer the files to be concatenated
#dx download -r --no-progress --lightweight *

find /home/dnanexus/src/results9/ -iname "*.filt" -print0 | xargs -0 awk '
NR == 1 {print "name_file\t" $0; next;}
FNR == 1 {next;}
{print FILENAME "\t" $0;}'
