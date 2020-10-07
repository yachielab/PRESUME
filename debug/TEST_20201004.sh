#!/bin/bash
# -*- cofing: utf-8 -*-
THIS_PATH=$(cd $(dirname $0); pwd)
PRESUME="${THIS_PATH}/../PRESUME.py"
OUTDIR="/Users/iwatano/Desktop/debug_PRESUME/after"

# run PRESUME
for i in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 5.0 10.0
do
presumeout="${OUTDIR}/presume_${i}"
mkdir $presumeout
python3 ${PRESUME} --output ${presumeout} -s $i --constant 0.1 -n 500 
done
echo "Done!"

