#!/bin/bash
#Defaults
############################################
raw=0
keepNormal=0

#Getopts
############################################
usageMessage="\nUsage: $0 [-h] -o outputFile [-r] [-k] inputfile1 ... inputfileN\n\nCommand line options:\n\t-h: help message\n\t-o: output file\n\t-r raw: It toggles using raw log2-ratios or corrected using cellularity and ploidy (default)\n\t-k keep: It toggles the elimination of normal (i.e. == ploidy) segments. Only works when -r is not active\n\n"
function usage {
    if [[ $# -eq 1 ]]
    then
        >&2 echo -e "\nERROR: $1"
        >&2 echo -e "$usageMessage"
    else
        >&2 echo -e "$usageMessage"
    fi

    exit 1
}


[ $# -eq 0 ] && usage
while getopts ":o:hkr" options
do
    case "${options}" in
        h)
            usage
            ;;
        o)
            outputFile=${OPTARG}
            ;;
		k)
			keepNormal=1
			;;
		r)
			raw=1
			;;
        *)
            usage "wrong input option"
            ;;
    esac
done
#After this, there would be positional arguments available
shift "$((OPTIND-1))"

if [[ ! -f "$lviProjectDNAScripts/configFile" ]]
then
	usage "The lviProjectDNAScripts environment variable has not been set or the configuration file configFile has not been set properly. Fix this issue before trying to execute this script again."
	exit 1
fi

source $lviProjectDNAScripts/configFile

echo -e "Generating the GISTIC2.0 output file $outputFile by parsing:"
echo -e "\t$@" | sed "s/ /\n\t/g"
echo -e "\n"

if [ "$raw" -eq 1 ]
then
	echo -e "Using raw log2ratios"
	[ "$keepNormal" -eq 1 ] && echo "WARNING: the -k option does not have an effect when combined with the -r option">&2
	
else
	if [ "$keepNormal" -eq 1 ]
	then
		echo -e "Using corrected log2ratios"
	else
		echo -e "Using corrected log2ratios and fixing Seg.CN to 0 in segments with round(absoluteCN) == round(ploidy)"
	fi
fi

mkdir -p $(dirname $outputFile)

echo -e "Sample\tChromosome\tStart Position\tEnd Position\tNum Markers\tSeg.CN" > $outputFile
for file in $@;
do
	sampleName=$(basename $file | sed "s/^\([^\.]*\.[^\.]*\).*/\\1/g")
	ploidy=$(basename $file | sed "s/.*\.p\([0-9]*\.*[0-9]*\).tsv/\\1/")
	if [ "$raw" -eq 0 ]
	then
		if [ "$keepNormal" -eq 1 ]
		then
			tail -n+2 $file  | awk -v sample="$sampleName" -v ploidy="$ploidy" 'BEGIN{OFS="\t"}{print(sample,$1,$2,$3,$4,log(($7<0?0.1:$7)/ploidy)/log(2))}' >> $outputFile
		else
			tail -n+2 $file  | perl -slane 'sub round {($n)=@_;int($n + $n/abs($n*2 || 1))};$rploidy=round($ploidy);$rACN=round($F[6]);if($rploidy==$rACN){$logR=0}else{$logR=log(($F[6]<0?0.1:$F[6])/$ploidy)/log(2)};print(join("\t",$sample,$F[0],$F[1],$F[2],$F[3],$logR))' -- -sample="$sampleName" -ploidy="$ploidy" >> $outputFile
		fi
	else
		tail -n+2 $file  | awk -v sample="$sampleName" 'BEGIN{OFS="\t"}{print(sample,$1,$2,$3,$4,$5)}' >> $outputFile
	fi
	
done
