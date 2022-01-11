#! /bin/bash

# all we have to do is exec python...
. /opt/conda/etc/profile.d/conda.sh
conda activate prfpy_analysis
echo "-------------- Activated conda environment ---------------"
exec python /run.py "$@"