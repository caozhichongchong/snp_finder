#!/bin/bash
set -e
keyCommand="ssh -i /Users/caozhichongchong/Desktop/AnniZhang/project/gut/c3ddb-cluster/linux/c3ddb-key"
echo rsync -e "$keyCommand" "$@" .
rsync -e "$keyCommand" "$@" .
