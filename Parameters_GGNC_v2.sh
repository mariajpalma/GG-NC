#!/bin/bash
export PYTHON_BIN=/export/apps/bioconda/bin/python3

show_help() {
    echo "Usage: bash Pipeline.sh -k <param1> -p <param2> -d <param3> -i <param4> -m <param5> -s <param6> -l <param7> -u <param8> -c <param9> -a <param10> -z <param11> -f <param12> -t <param13> -r <param14> -T <param15> [-L]"
    echo "Options:"
    echo "  -k <param1>    Kind: IBD, PCA or GRM"
    echo "  -p <param2>    Path of your files"
    echo "  -d <param3>    Name of your data files"
    echo "  -i <param4>    Name of your info file"
    echo "  -m <param5>    Maximum value"
    echo "  -s <param6>    Steps in the log10 space to explore"
    echo "  -l <param7>    Lambda value to explore"
    echo "  -u <param8>    Prune option"
    echo "  -c <param9>    Minimun individuals in a community"
    echo "  -a <param10>   Lower limit of the log10 space to explore"
    echo "  -z <param11>   Upper limit of the log10 space to explore"
    echo "  -f <param12>   Skip network plots"
    echo "  -t <param13>   Number of cores"
    echo "  -r <param14>   (Optional) Set randome seed"
    echo "  -T <param15>   (Optional) File with pairwise mean TMRCA"
    echo "  -L             Use the Leiden algorithm instead of Louvain"
    exit 1
}


# Initialize variables
use_leiden=false
random_seed=""
TMRCA_file=""

#Reading parameters
while getopts ":k:p:d:i:m:s:l:u:c:a:z:f:t:r:T:L" option; do
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
    a)
      lower_limit="$OPTARG"
      ;;
    z)
      upper_limit="$OPTARG"
      ;;
    f)
      figures="$OPTARG"
      ;;
    t)
      n_cores="$OPTARG"
      ;;
    r)
      random_seed="$OPTARG"
      ;;
    T)
      TMRCA_file="$OPTARG"
      ;;
    L)
      use_leiden=true
      ;;
    *)
      show_help
      ;;
  esac
done



if [ -z "$kind" ] || [ -z "$path" ] || [ -z "$data_files" ] || [ -z "$info_file" ] || [ -z "$max" ] || [ -z "$steps" ] || [ -z "$lambda" ] || [ -z "$prune" ] || [ -z "$min_comms" ] || [ -z "$upper_limit" ] || [ -z "$lower_limit" ] || [ -z "$figures" ] || [ -z "$n_cores"]; then
    echo "Some parameters are missing."
    show_help
fi


if [ -n "$TMRCA_file" ]; then
    TMRCA_path="${path%/}/$TMRCA_file"

    if [ ! -f "$TMRCA_path" ]; then
        echo "The TMRCA file does not exist: $TMRCA_path"
        exit 1
    fi
fi

DATA_path="${path%/}/$data_files"

if [ ! -f "$DATA_path" ]; then
    echo "The data file does not exist: $DATA_path"
    exit 1
fi

INFO_path="${path%/}/$info_file"

if [ ! -f "$INFO_path" ]; then
    echo "The info file does not exist: $INFO_path"
    exit 1
fi


echo "The genetic metric is: $kind"
echo "The path is: $path"
echo "The input file name is: $data_files"
echo "The info file name is: $info_file"
echo "The max value is: $max"
echo "The number of steps is: $steps"
echo "The lambda value is: $lambda"
echo "The option for prune is: $prune"
echo "The minnimun number of detected communities is: $min_comms"
echo "The lower limit of the log10 space to explore is: $lower_limit"
echo "The upper limit of the log10 space to explore is: $upper_limit"
echo "The generation of network plots is set to: $figures"
echo "The number of cores is set to: $n_cores"

if [ -n "$random_seed" ]; then
    echo "The random seed was set to: $random_seed"
else
    echo "No random seed was provided."
fi

if [ -n "$TMRCA_file" ]; then
    echo "TMRCA file passed with -T: $TMRCA_file"
fi

if [ "$use_leiden" = true ]; then
    echo "Using Leiden algorithm."
else
    echo "Using Louvain algorithm."
fi


# Build arg array

args=(
  "$kind"
  "$path"
  "$data_files"
  "$info_file"
  "$max"
  "$steps"
  "$lambda"
  "$prune"
  "$min_comms"
  "$lower_limit"
  "$upper_limit"
  "$figures"
  "$n_cores"
)


if [ -n "$random_seed" ]; then
    args+=("$random_seed")
fi

if [ -n "$TMRCA_file" ]; then
    args+=("$TMRCA_file")
fi

if [ "$use_leiden" = true ]; then
    args+=("-L")
fi

# Call Rscript
Rscript CommunityNetwork-Analysis_v07052025.R "${args[@]}"
