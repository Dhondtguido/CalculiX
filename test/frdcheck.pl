#!/usr/bin/perl

#
#   usage: frdcheck.pl file
#          compares file.frd with file.frd.ref
#

$file=$ARGV[0];

#
#   determining the maximum in each data block
#

open FRDREF, "<$file.frd.ref";
$k=0;
$xmax[0]=0;
$active=0;
while(<FRDREF>){
    if(/^ -5/){
	$k++;
	$xmax[$k]=0.;
	$active=1;
	next;
    }

    if(/^ -3/){
	$active=0;
    }

    if($active==1){
	for($i=1;$i<=$#fieldsref;$i++){$fieldsref[$i]=0;}
	$count=tr/ .eE+\-0-9//;
	for($i=1;$i<=($count-13)/12;$i++){
	    $fieldsref[$i]=substr($_,1+12*$i,12);
	}
	$n=($count-13)/12;
	for($i=1;$i<=$n;$i++){
	    if($fieldsref[$i] =~ /\./){
		if(abs($fieldsref[$i])>$xmax[$k]){
		    $xmax[$k]=abs($fieldsref[$i]);
		}
	    }
	}
    }
}
close FRDREF;

open FRD, "<$file.frd";
open FRDREF, "<$file.frd.ref";

$k=0;
$j=0;
$maxerror2=0.;
$active=0;

while(<FRDREF>){
#
#   lines starting with -5 are not compared
#
    if(/^ -5/){
	$k++;
	$active=0;
    }
#
#   comparison stops at a line starting with -3
#
    if(/^ -3/){
	$active=0;
    }
#
#   incrementing the line number
#
    $j++;
#
#   lines in active blocks in the reference file are split
#
    if($active==1){
	$count=tr/ .eE+\-0-9//;
	for($i=1;$i<=($count-13)/12;$i++){
	    $fieldsref[$i]=substr($_,1+12*$i,12);
	}
	$n=($count-13)/12;
    }

    while(<FRD>){
#
#       jump for inactive lines
#
	if($active==0){last;}
#
#       lines in active blocks in the frd file are split
#
	$count=tr/ .eE+\-0-9//;
	for($i=1;$i<=($count-13)/12;$i++){
	    $fields[$i]=substr($_,1+12*$i,12);
	}
#
#       comparing the values in the same line
#
	for($i=1;$i<=$n;$i++){
	    if(abs($fieldsref[$i])<1.e-10){next;}
	    if($xmax[$k]>1.e-10){
		$error2=abs(($fieldsref[$i]-$fields[$i])/$xmax[$k]);
	    } 
	    else {
		$error2=0.;
	    }
#
#           error2 is the relative error w.r.t. the maximum reference
#           value in the same block (displacements, stresses...)
#
	    if(($error2>1.e-3)&&($error2>$maxerror2)){
		$value2=$fields[$i];
		$valueref2=$fieldsref[$i];
		$line2=$j;
		$maxerror2=$error2;
		$largestvalue=$xmax[$k];
	    }
	}
	last;
    }
#
#   the line after a line starting with -5 may be compared
#
    if(/^ -5/){
	$active=1;
    }

}

if($maxerror2>0.){
   $maxerror2=100.*$maxerror2;
    printf "deviation in file $file.frd\n";
    printf "line: %d reference value: %e value: %e \n",$line2,$valueref2,$value2;
    printf "            absolute error: %e \n",abs($value2-$valueref2);    
    printf "            largest value within same block: %e \n",$largestvalue;
    printf "            relative error w.r.t. largest value within same block: %f \%\n\n",$maxerror2;
}

close FRD;
close FRDREF;

	
