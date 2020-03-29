#!/bin/bash
# -*- cofing: utf-8 -*-
THIS_PATH=$(cd $(dirname $0); pwd)
PRESUME="${THIS_PATH}/../PRESUME.py"
OUTDIR=${THIS_PATH}/test
TREE="${THIS_PATH}/../example/example_1/PRESUMEout/PRESUMEout.nwk"
# run PRESUME
mkdir ${OUTDIR}
python3 ${PRESUME} --output ${OUTDIR} --tree  ${TREE}
cat ${TREE} | sed -e 's/[0-9]//g' > ${OUTDIR}/reference1
cat ${OUTDIR}/PRESUMEout_from_tree/test.nwk | sed -e 's/[0-9]//g' > ${OUTDIR}/reference2
diff -s ${OUTDIR}/reference1 ${OUTDIR}/reference2
echo "Done!"