#!python
import os, sys

def usage():
  print("usage: sam-select.py <alignments.sam> <contig name> <position>", file=sys.stderr)
  print("outputs alignments that intersect the given contig name and position", file=sys.stderr)
  sys.exit(1)

def run(alignmentsFile, filterContigName, filterPosition):
  for line in open(alignmentsFile):
    if line.startswith("@"):
      continue # comment
    columns = line.split("\t")
    queryContigName = columns[2]
    queryStart = int(columns[3])
    queryText = columns[9]
    queryLength = len(queryText)
    queryEnd = queryStart + queryLength # this isn't quite right in case of an indel but it is probably close enough
    if filterContigName == queryContigName:
      if filterPosition >= queryStart and filterPosition <= queryEnd - 1:
        print(line, end='')

def main(args):
  if len(args) < 4:
    usage()
  alignmentsFile = args[1]
  filterContigName = args[2]
  filterPosition = int(args[3])
  run(alignmentsFile, filterContigName, filterPosition)

if __name__ == "__main__":
  main(sys.argv)
