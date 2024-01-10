#!/usr/bin/perl

#
#   usage: datcheck.pl file
#          compares file.dat with file.dat.ref
#

$file=$ARGV[0];

#
#   determining the maximum in each data block
#

open DATREF, "<$file.dat.ref";
$k=0;
$xmax[0]=0;
while(<DATREF>){
    if(/^$/){
	$k++;
	$xmax[$k]=0.;
	next;
    }
    @fieldsref=split /\s+/, $_;
    $n=$#fieldsref;
    for($i=1;$i<=$n;$i++){
#    for($i=3;$i<=$n;$i++){
	if($fieldsref[$i] =~ /\./){
	    if(abs($fieldsref[$i])>$xmax[$k]){
		$xmax[$k]=abs($fieldsref[$i]);
	    }
	}
    }
#    print $k,$xmax[$k],"\n";
}
close DATREF;

open DAT, "<$file.dat";
open DATREF, "<$file.dat.ref";

$k=0;
$j=0;
$maxerror1=0.;
$maxerror2=0.;

# checks for a participation factor block in a *frequency calculation :
# should not be compared (dependent on the cyclic sector)

$participation=0;

while(<DATREF>){
    if(/^$/){
	$k++;
    }
    if(/P A R T I C I P A T I O N   F A C T O R S   F O R   F R E Q U E N C Y/){
    }elsif(/P A R T I C I P A T I O N   F A C T O R S/){
	$participation=1;
    }elsif(/T O T A L   E F F E C T I V E   M A S S/){
	$participation=0;
    }
    $j++;
    @fieldsref=split /\s+/, $_;
    while(<DAT>){
	@fields=split /\s+/, $_;
	$n=$#fieldsref;
	if(/time/){last;}
	if($participation==1){last;}
	for($i=1;$i<=$n;$i++){
	    if(abs($fieldsref[$i])<1.e-10){next;}
	    $error1=abs(($fieldsref[$i]-$fields[$i])/$fieldsref[$i]);
	    if($xmax[$k]>1.e-10){
		$error2=abs(($fieldsref[$i]-$fields[$i])/$xmax[$k]);
	    } 
	    else {
		$error2=0.;
	    }
#
#           error1 is the relative error w.r.t. the reference value
#           at the same position
#
	    if(($error1>1.e-3)&&($error1>$maxerror1)&&(abs($fieldsref[$i])>1.e-5)){
		$value1=$fields[$i];
		$valueref1=$fieldsref[$i];
		$line1=$j;
		$maxerror1=$error1;
	    }
#
#           error2 is the relative error w.r.t. the maximum reference
#           value in the same block (displacements, stresses...)
#
	    if(($error2>1.e-3)&&($error2>$maxerror2)&&
               (($fields[$i]>1.e-8)||($xmax[$k]>1.e-8))){
		$value2=$fields[$i];
		$valueref2=$fieldsref[$i];
		$line2=$j;
		$maxerror2=$error2;
		$largestvalue=$xmax[$k];
	    }
	}
	last;
    }
}

#if($maxerror1>0.){
#    $maxerror1=100.*$maxerror1;
#    print "deviation in file $file.dat\n";
#    print "line: $line1 reference value: $valueref1 value: $value1 \n";
#    printf "            absolute error: %f \n",abs($value1-$valueref1);    
#    printf "            relative error w.r.t. reference value at same position: %f \%\n\n",$maxerror1;
#}

if($maxerror2>0.){
   $maxerror2=100.*$maxerror2;
    print "deviation in file $file.dat\n";
    printf "line: %d reference value: %e value: %e \n",$line2,$valueref2,$value2;
    printf "            absolute error: %e \n",abs($value2-$valueref2);    
    printf "            largest value within same block: %e \n",$largestvalue;
    printf "            relative error w.r.t. largest value within same block: %f \%\n\n",$maxerror2;
}

close DAT;
close DATREF;

	
