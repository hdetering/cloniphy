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

#SBATCH -J cloniphy              # Job name
#SBATCH -o cloniphy.%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p cola-corta    #-p shared --qos=shared
#SBATCH -n 10                    # Total number of tasks
#SBATCH -t 10:00:00              # Run time (hh:mm:ss)

CLONIPHY=$HOME/code/cloniphy/build/src/demo
#CONFIG=$HOME/code/cloniphy/config.min.cesga.yml
CONFIG=$LUSTRE/cloniphy/config/config.ngsws.yml
OUTDIR=$LUSTRE/cloniphy/output
THREADS=$SLURM_NTASKS

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
echo "[INFO] $CLONIPHY -c $CONFIG -o $TMPDIR/$PRJ 1>$TMPDIR/cloniphy.out 2>$TMPDIR/cloniphy.err"
time $CLONIPHY -c $CONFIG -o $TMPDIR/$PRJ 1>$TMPDIR/cloniphy.out 2>$TMPDIR/cloniphy.err

# Post-process output.
echo "[INFO] Post-processing output..."

echo "[INFO]   cp $CONFIG $TMPDIR/$PRJ"
cp $CONFIG $TMPDIR/$PRJ

echo "[INFO]   mv $TMPDIR/cloniphy.{out,err} $TMPDIR/$PRJ"
mv $TMPDIR/cloniphy.{out,err} $TMPDIR/$PRJ

echo "[INFO]   samtools faidx $TMPDIR/$PRJ/ref.fa"
samtools faidx $TMPDIR/$PRJ/ref.fa

for f in $TMPDIR/$PRJ/bam/sample*.sam; do
  bam="${f%sam}bam"
  echo "[INFO]  samtools sort -@ $THREADS $f > $bam && samtools index $bam"
  samtools sort -@ $THREADS $f > $bam && samtools index $bam
  if [ -f "$bam" ] && [ -s "$bam" ]; then
    echo "         rm $f"
    rm $f
  fi 
done

# Move output to target dir.
echo
echo "cp -r $TMPDIR/$PRJ $OUTDIR/"
cp -r $TMPDIR/$PRJ $OUTDIR/
