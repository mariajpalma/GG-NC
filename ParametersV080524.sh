#!/bin/bash

show_help() {
    echo "Usage: bash Pipeline.sh -k <param1> -p <param2> -d <param3> -i <param4> -m <param5> -s <param6> -l <param7> -u <param8> -c <param9>"
    echo "Options:"
    echo "  -k <param1>    Kind: IBD, PCA or GRM"
    echo "  -p <param2>    Path of your files"
    echo "  -d <param3>    Name of your data files"
    echo "  -i <param4>    Name of your info file"
    echo "  -m <param5>    Maximum value"
    echo "  -s <param6>    Steps"
    echo "  -l <param7>    Lambda value to explore"
    echo "  -u <param8>    Prune option"
    echo "  -c <param9>    Minimun individuals in a community"
    exit 1
}

#Reading parameters
while getopts ":k:p:d:i:m:s:l:u:c:" option; do
  case $option in
    k)
      kind="$OPTARG"
      ;;
    p)
      path="$OPTARG"
      ;;
    d)
      data_files="$OPTARG"
      ;;
    i)
      info_file="$OPTARG"
      ;;
    m)
      max="$OPTARG"
      ;;
    s)
      steps="$OPTARG"
      ;;
    l)
      lambda="$OPTARG"
      ;;
    u)
      prune="$OPTARG"
      ;;
    c)
      min_comms="$OPTARG"
      ;;
    *)
      show_help
      ;;
  esac
done

if [ -z "$kind" ] || [ -z "$path" ] || [ -z "$data_files" ] || [ -z "$info_file" ] || [ -z "$max" ] || [ -z "$steps" ] || [ -z "$lambda" ] || [ -z "$prune" ] || [ -z "$min_comms"]; then
    echo "Some parameters are missing."
    show_help
fi

echo "The genetic metric is: $kind"
echo "The path is: $path"
echo "The input file name is: $data_files"
echo "The info file name is: $info_file"
echo "The max value is: $max"
echo "The number of steps is: $steps"
echo "The lambda value is: $lambda"
echo "The option for prune is: $prune"
echi "The minnimun number of detected communities is: $min_comms"

### Call Rscript 
#conda activate r422

Rscript CommuntyNetwork-Analysis.R $kind $path $data_files $info_file $max $steps $lambda $prune $min_comms
