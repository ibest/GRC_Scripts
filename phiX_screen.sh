#!/bin/bash

# A POSIX variable
OPTIND=1  # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
threads=4
pe1=NA
pe2=NA
se="-"


while getopts "h?t:1:2:U:" opt; do
    case "$opt" in
    h|\?)
        echo "there is no help"
        exit 0
        ;;
    v)  threads=$OPTARG
        ;;
    1)  pe1=$OPTARG
        ;;
    2)  pe2=$OPTARG
        ;;
    U)  se=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

btargs=$@

bt_command="bowtie2 -I 0 -X 1500 --very-sensitive-local -p $threads -x phiX.index -q "
if [ $se = "-" ] ; then # expect paired end reads
   if [[ ${pe1: -3} == ".gz" ]] ; then # set up named pipes
      mkfifo pe1 pe2
      declare -a prefix=("zcat $pe1 &> pe1 &" "zcat $pe2 &> pe2 &")
      suffix=" -1 pe1 -2 pe2"
   else
      suffix=" -1 $pe1 -2 $pe2"
   fi
else
   if [[ ${se: -3} == ".gz" ]] ; then # set up named pipes
      mkfifo se 
      declare -a prefix=(" zcat $se > se &")
      suffix=" -U se"
   else
      suffix=" -U $se"
   fi
fi
#if [[ -z "${prefix[0]}" ]]; then
#   for cmd in "${prefix[@]}" ; do
#      echo "$cmd"
#      $cmd
#   done
#fi
echo "${prefix[0]}"
${prefix[0]}
echo "${prefix[1]}"
${prefix[1]}

source /usr/modules/init/bash
module load bowtie2 grc

echo "The command:"
echo "$bt_command $suffix | samtools view -bS "
$bt_command $suffix
# | samtools view -bS 

if [ -f pe1 ] ; then rm pe1; fi
if [ -f pe2 ] ; then rm pe2; fi


