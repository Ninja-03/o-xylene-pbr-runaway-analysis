#!/bin/sh
cd /mnt/user-data/outputs/extension_study
for m in ds2 delasa cc; do python3 pi121_run.py $m > pi121_$m.log 2>&1; done
echo ALLDONE > pi121_done.flag
