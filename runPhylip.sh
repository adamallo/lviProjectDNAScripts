#!/bin/bash

#Run seqboot
for i in P[0-9][0-9].phy; do cp $i infile; /Applications/phylip/seqboot.app/Contents/MacOS/seqboot < ../../scripts/seqbootCMD > $(echo $i | sed "s/.phy/_seqboot.log/g"); cp outfile $(echo $i | sed "s/.phy/_bootstrap.phy/");done

#Run dollop for all bootstrapped trees
for i in *bootstrap.phy; do cp $i infile; /Applications/phylip/dollop.app/Contents/MacOS/dollop < ../../scripts/dollopCMDBoost > $(echo $i | sed "s/.phy/_dollop.log/g"); cp outtree $(echo $i | sed "s/.phy/_boostrapDollop.trees/");done

#Run dollop for the full datasets
for i in P[0-9][0-9].phy; do cp $i infile; /Applications/phylip/dollop.app/Contents/MacOS/dollop < ../../scripts/dollopCMD > $(echo $i | sed "s/.phy/_dollopFull.log/g"); cp outtree $(echo $i | sed "s/.phy/_dollopFull.trees/");done

#Transform phylip bootstrap outputs with weights to trees
for file in *_bootstrap_boostrapDollop.trees;do reweight=$(grep -o -e '\[.*\]' $file | sort | uniq | head -n 1 | perl -ne 's/\[(.*)\]/$1/;printf "%.0f",1/$_'); cat $file | perl -pe '/,\n/ and chomp' | perl -ne 'BEGIN{our $maxNlines=6}{if($_ =~ m/\[/){(my $temp = $_) =~s/^.*\[(.*)\].*$/$1/g;$nlines=$maxNlines/int(sprintf("%.0f",1/$temp));($line=$_)=~s/\[.*\]//}else{$line=$_;$nlines=$maxNlines};for(my $i=0;$i<$nlines;$i++){print $line}}' > $(echo $file | sed "s/_bootstrap_boostrapDollop/_bootstrap_boostrapDollopFixed/g");done

