#!/bin/bash

# A POSIX variable
OPTIND=1  # Reset in case getopts has been used previously in the shell.

# Initialize variables:
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
    o)  output=$OPTARG
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

#check for the bowtie2 index files
if [ -f /data/phiX/phiX.index ] ; then 
   bowtie_index=/data/phiX/phiX.index
else
   wget -O phiX.tmp.fasta "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=9626372"
   if [ $? -eq 0 ] ; then
     bowtie2-build phiX.tmp.fasta phiX.bowtie2_tmp_index
     bowtie_index=phiX.bowtie2_tmp_index
   else
     echo "ERROR: No copy of phiX sequence available, and download from NCBI failed. Goodbye"
     exit 1
   fi
fi

#build the command string
bt_command="bowtie2 -I 0 -X 1500 --very-sensitive-local -p $threads -x $bowtie_index -q "
if [ $se = "-" ] ; then # expect paired end reads
      suffix=" -1 $pe1 -2 $pe2"
else
      suffix=" -U $se"
fi
source /usr/modules/init/bash
module load bowtie2 grc/2.0

#echo "The command:"
#echo "$bt_command $suffix | samtools view -bS "
#$bt_command $suffix | samtools view -bS - > btout.bam
$bt_command $suffix | extract_unmapped_reads.py -o $output_noContaminant

# clean up
if [ -f phiX.bowtie2_tmp_index ] ; then rm phiX.tmp.fasta phiX.bowtie2_tmp_index.* phiX.bowtie2_tmp_index; fi

exit 0
