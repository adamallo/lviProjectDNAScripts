#!/bin/bash

if [[ ! -f "$lviProjectDNAScripts/configFile" ]]
then
    echo "The lviProjectDNAScripts environment variable has not been set or the configuration file configFile has not been set properly. Fix this issue before trying to execute this script again."
    exit 1
fi

source $lviProjectDNAScripts/configFile

#Run seqboot
for i in *.phy; do cp $i infile; $phylipDir/seqboot < $scriptsDir/seqbootCMD > $(echo $i | sed "s/.phy/_seqboot.log/g"); cp outfile $(echo $i | sed "s/.phy/_bootstrap.phy/");done

#Run dollop for all bootstrapped sequences
for i in *bootstrap.phy; do cp $i infile; $phylipDir/dollop < $scriptsDir/dollopCMDBoost > $(echo $i | sed "s/.phy/_dollop.log/g"); cp outtree $(echo $i | sed "s/.phy/_bootstrapDollop.trees/");done

#Run dollop for the full datasets
for i in $(ls *.phy | grep -v bootstrap); do cp $i infile; $phylipDir/dollop < $scriptsDir/dollopCMD > $(echo $i | sed "s/.phy/_dollopFull.log/g"); cp outtree $(echo $i | sed "s/.phy/_dollopFull.trees/");done

#Transform phylip bootstrap outputs with weights to trees
for file in *_bootstrap_bootstrapDollop.trees;do reweight=$(grep -o -e '\[.*\]' $file | sort | uniq | head -n 1 | perl -ne 's/\[(.*)\]/$1/;printf "%.0f",1/$_'); cat $file | perl -pe '/,\n/ and chomp' | perl -ne 'BEGIN{our $maxNlines=6}{if($_ =~ m/\[/){(my $temp = $_) =~s/^.*\[(.*)\].*$/$1/g;$nlines=$maxNlines/int(sprintf("%.0f",1/$temp));($line=$_)=~s/\[.*\]//}else{$line=$_;$nlines=$maxNlines};for(my $i=0;$i<$nlines;$i++){print $line}}' > $(echo $file | sed "s/_bootstrap_bootstrapDollop/_bootstrap_bootstrapDollopFixed/g");done

