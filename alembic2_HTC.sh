#!/bin/bash

# prepare a working directory (named as outDir since the output files will be stored here)

# outDir = ${1%.*_*.*.*.*.*.*}

outDir=$1"_dir"
mkdir $outDir

mv $1 $outDir
mv calculate_gap_lengths.pl $outDir

# goto the working directory ($outDir) and run the perl script
cd $outDir
cat $1 | ./calculate_gap_lengths.pl > $1".zanfona.joins.csv"

# remove those files which are no longer needed; i.e. don't need to transfer them back to submission node
rm $1
rm calculate_gap_lengths.pl 

# the rest of the files inside $outDir are those needed to be transferred back!
# so, in your condor sumission file, use the "transfer_output_files = ... "

