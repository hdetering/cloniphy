#!/bin/bash
#==============================================================================
# title           :run_cloniphy.cesga.sh
# description     :Runs CloniPhy in an interactive session.
# author          :Harald Detering (harald.detering@gmail.com
# date            :2017-10-25
# version         :0.1
# usage	          :bash run_cloniphy.cesga.sh <project_name>
# notes           :Uses computing node's local disk for faster file I/O.
#                  Results will be moved to target dir specified in the script.
#==============================================================================

CLONIPHY=$HOME/code/cloniphy/build/src/demo
#CONFIG=$HOME/code/cloniphy/config.min.yml
CONFIG=$LUSTRE/cloniphy/config/config.ngsws.yml
OUTDIR=$LUSTRE/cloniphy/output

if [ $# -lt 1 ]; then
  (>&2 echo)
  (>&2 echo "usage: $0 <project_name>")
  (>&2 echo)
  exit 1
fi


# load required modules
module load gcc/5.3.0
module load boost/1.61.0

# This will be the name of the output folder.
PRJ="$1"

# Run CloniPhy, writing output to local disk.
echo "[INFO] Running CloniPhy..."
time $CLONIPHY -c $CONFIG -o $TMPDIR/$PRJ 1>$TMPDIR/cloniphy.out 2>$TMPDIR/cloniphy.err

# Post-process output.
echo "[INFO] Post-processing output..."
for f in $TMPDIR/$PRJ/bam/sample*.sam; do
  echo "[INFO]   $f"
  bam="${f%sam}bam"
  samtools sort $f > $bam && samtools index $bam
  if [ -f "$bam" ] && [ -s "$bam" ]; then
    rm $f
  fi 
done

# Move output to target dir.
mv $TMPDIR/$PRJ $OUTDIR/
