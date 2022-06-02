set -e

function usage() {
  echo "usage: $0 <filepath> <contigName> <startIndex> <length>"
  exit 1
}


filepath="$1"
contigName="$2"
startIndex="$3"
length="$4"

if  [ "$filepath" == "" ]; then
  usage
fi

if [ "$contigName" == "" ]; then
  usage
fi

if [ "$startIndex" == "" ]; then
  usage
fi

if [ "$length" == "" ]; then
  usage
fi

endIndex="$((startIndex + length - 1))"

#for i in $(seq 1 $startIndex); do
#  echo -n " "
#done
grep "$contigName" -A 1 "$filepath" | tail -n 1 | head -c "$endIndex" | tail -c "$length"
echo
