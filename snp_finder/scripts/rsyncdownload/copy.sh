#!/bin/bash

fromDir=anniz44@c3ddb01.mit.edu:/scratch/users/anniz44/Metagenomes/BN10_MG/assembly
keyCommand="ssh -i /Users/caozhichongchong/Desktop/AnniZhang/project/gut/c3ddb-cluster/linux/c3ddb-key"

echo finding files
rsync -e "$keyCommand"  $fromDir/*.gz . --dry-run --stats -avzm > log

echo parsing list of files
cat log | grep gz | sed "s|^|$fromDir/|" > latest.gz-list

echo calculating files per worker
numFiles="$(cat latest.gz-list | wc -l)"
numWorkers="$numFiles"
numFilesPerWorker="$(($numFiles / $numWorkers))"

echo building rsync command
cat latest.gz-list | xargs -n "$numFilesPerWorker" echo ./rsync-dot.sh  > commands.sh

echo writing scripts
rm -f x*
split -l 1 commands.sh
chmod u+x x*

echo Now run "x*" in parallel
