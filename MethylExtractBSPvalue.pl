#!/usr/bin/perl -w

=about
###########################################################################
###########################################################################

  *****************************************  
  *****  MethylExtract (version 1.9.0) ****  
  *****************************************  


 Computational Epigenomics and Bioinformatics
  Dept. of Genetics & Inst. of Biotechnology 
             University of Granada           

         Web: http://bioinfo2.ugr.es/        

 This program is Copyright (C) 2012-16:
 Ricardo Lebrón (rlebron@ugr.es), Guillermo Barturen (bartg01@gmail.com), Antonio Rueda (aruemar@gmail.com), José L. Oliver (oliver@ugr.es), Michael Hackenberg (hackenberg@ugr.es)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 For questions, feedback, etc. please contact to:
 Ricardo Lebrón (rlebron@ugr.es), Guillermo Barturen (bartg01@gmail.com), José L. Oliver (oliver@ugr.es, Michael Hackenberg (mlhack@gmail.com)

 To see the options, please launch MethylExtract without any command line arguments


############################################################################
############################################################################
=cut

use strict;

#################
#Default options#
#################
my %optionsDefault = (
	inFile => "NA", outFile => "NA",
 	BSCR => "NA",
 	errorInterval => "0.2", FDR => "NA",
);
#################

#Getting options
my ($infile,$outfile,$prob,$interval,$fdr);
($infile,$outfile,$prob,$interval,$fdr)=
&GetOptions($prob,$interval,$fdr);

print "\n################# Running MethylExtractBSPvalue v1.9 #################\n\n";

open(I,$infile) or die "could not open $infile";
open(O,">$outfile") or die "could not open $outfile";

my $cl = 0;
my $count=0;
my $cytosineCalculated=0;
my $countCalculated=0;
print "Calculating p-values for each cytosine methylation value\n";
while(my $z = <I>){
	if (substr($z,0,1) eq "#") {
		chomp($z);
		print O "$z\tBS Error Probability\n";
		next;
	}
	$count++;
	my @f = split("\t",$z);		
	if(@f > $cl){
		$cl = scalar @f;
	}
	my $trials=0;
	my $mC=0;
	if ($f[3] ne ".") {
		$trials+=$f[4];
		$mC+=$f[3];
	}
	if ($f[6] ne ".") {
		$trials+=$f[7];
		$mC+=$f[6];
	}
	#my $observedLevel = $mC/$trials;
	chomp($z);
	if ($mC>$trials) {
		print "Warning: methylated cytosines are higher than coverage at $z, the value has been discarded!!!\n";
		next;
	}
	my $maxAllowedBSfail = &allowedBSfailures($trials, $mC, $interval);		
	my $tempP=0;
	# CDF
	my $p=0;
	if ($mC>0) {
		for(my $i=0;$i<$maxAllowedBSfail;$i++){
			my $probTemp= &Binomial($mC, $i, $prob);
			$tempP +=  $probTemp;
		}
		$p=1-$tempP;
		if ($p<0) {$p=0}
	}
	else {}
	$p=sprintf("%.3e",$p);
	print O "$z\t$p\n";
	$cytosineCalculated++;
	if ($cytosineCalculated==100000) {
		$countCalculated+=$cytosineCalculated;
		$cytosineCalculated=0;
		print "$countCalculated calculated so far...\n";
	}
}
close(O);
	
if($fdr ne "NA"){
	print "\nCalculating FDR\n";
	$cl=$cl+1;
	print "will order for column $cl \n";
	my $code = system("sort -k $cl -g -o $outfile"."s ".$outfile);
	my $kmax = &fdr($outfile."s",$count,$cl-1,$outfile,$outfile.".noSig");
	unlink($outfile."s");
	printf ("found %d out of %d significant \n",$kmax,$count);
}
	
#############SUBPROCESS###################################################

sub GetOptions{
my %opts;
#Checking General Arguments & Help
if (@ARGV) {
	foreach (@ARGV) {
		my @GetOpt=split(/=/,$_);
		$opts{$GetOpt[0]}=$GetOpt[1];
	}
	#Checking input options
	foreach my $keyOpts (keys %opts){
		if (!$optionsDefault{$keyOpts}) {die "$keyOpts is not an accepted parameter, please check spelling and case sensitive\n";}
		else {}
	}
	#Checking input file
	if($opts{inFile}){
		$infile=$opts{inFile};
		#Checking input file
		my $checkAlign="N";
		if (-e $infile) {	
			$infile =~ m/(\w+)$/;
			if ($1 eq "output") {$checkAlign="Y";}
			else {}
			if ($checkAlign eq "N") {die "The input file doesn't seem to be a MethylExtract output file\n";}
			else {}
		}
		else {die "Cannot find input file: $infile\n";}
	}
	else {die "Use inFile=[input file] to specify the input file\n";}
	#Checking output file
	if ($opts{outFile}) {
		$outfile=$opts{outFile};
	}
	else {$outfile="$infile.prob";}
	#Checking error probability
	if (defined($opts{BSCR})) {
		$_[0]=$opts{BSCR};
		if ($_[0]=~m/[0-9]+/) {
			$_[0]=~s/,/\./;
			if (int($_[0])>1) {die "Bisulfite conversion rate must be between 0-1\n";}
			else {$_[0]=1-$_[0]}
		}
		else {die "Bisulfite conversion rate must be between 0-1\n";}
	}
	else {die "Use BSCR=[Bisulfite conversion rate] \n";}
	#Checking error interval
	if (defined($opts{errorInterval})) {
		$_[1]=$opts{errorInterval};
		if ($_[1]=~m/[0-9]+/) {
			$_[1]=~s/,/\./;
			if (int($_[1])>1) {die "Methylation interval must be between 0-1\n";}
			else {}
		}
		else {die "Methylation interval must be between 0-1\n";}
	}
	else {$_[1]=0.2;}
	#Checking FDR input
	if (defined($opts{FDR})) {
		$_[2]=$opts{FDR};
		if ($_[2]=~m/[0-9]+/) {
			$_[2]=~s/,/\./;
			if (int($_[2])>1) {die "FDR value must be between 0-1\n";}
			else {}
		}
		else {die "FDR value must be between 0-1\n";}
	}
	else {$_[2]="NA";}
	
}
else {
	print "\n################   MethylExtractErrorProbability   ###############\n";
	print "######################   Command-line help   #####################\n\n";
	print "Launch as:\n  perl MethylExtractBSPvalue.pl inFile=<input file> BSCR=<Bisulfite conversion rate> [OPTIONS]\n\n";
	print "Optional parameters:\n";
	print "  outFile=<Output file> [default: inFile.prob]\n";
	print "  errorInterval=<Error interval allowed> [default: 0.2]\n";
	print "  FDR=<False discovery rate allowed> [default: NA]\n";
	die "\n";
}
return($infile,$outfile,$_[0],$_[1],$_[2]);
}

## get the last p-value that is statistically significant	
## $_[0] --> the sorted file
## $_[1] --> the number of tests
## $_[2] --> the column with the p-value
## $_[3] --> the outfile with significant positions
## $_[4] --> not significant positions
sub fdr{
	open(I,"$_[0]") or die "could not open $_[0] within getK";
	open(SIG,">$_[3]") or die "could not open $_[3]";
	open(NOSIG,">$_[4]") or die "could not open $_[4]";
	my $k = 1;
	while(my $z = <I>) {
		if (substr($z,0,1) eq "#") {
			print SIG "$z";
			print NOSIG "$z";
			next;
		}
		my @f = split("\t",$z);		
		if($fdr*$k/$_[1] >= $f[$_[2]]){
			print SIG $z;
			$k++;
		}
		else{
			print NOSIG $z;
		}
	}
	close(SIG);
	close(NOSIG);
	close(O);
	return $k-1;
}

### 0--> trials(read coverage);; 1-->mC;; 2-->interval
sub allowedBSfailures{
	my $back = 0;
	my $observed = $_[1]/$_[0];
	for(my $i = 0; $i <= $_[1]; $i++){
		my $act = ($_[1]-$i)/$_[0];
		if(($observed - $act) le $_[2]) {
			$back = $i;
		}
		else{
			return $back;
		}
	}
	return $back;
}

## 0--> trials, 1--> number of successes, 2--> probability
sub Binomial{
	my $temp = &BinomialCoeff($_[0],$_[1]) + $_[1]*log($_[2]) + ($_[0]-$_[1])*log(1.0-$_[2]);
	return exp($temp);
}

## 0->over; 1-->under
sub BinomialCoeff{
	my $over = $_[0];
	my $under = $_[1];
	if($under > $over){
		return -1;
	}
	elsif($under == 0){
		return 0;
	}
	elsif($under < $over - $under){
		# first sum
		my $first = 0;
		my $second = 0;
		for(my $i = $over-$under+1; $i <= $over; $i++){
			$first += log($i);
		}
		# $second sum
		for(my $i = 1; $i <= $under; $i++){
			$second += log($i);
		}
		return $first-$second;
	}	
	else{
		my $first = 0;
		my $second = 0;
		for(my $i = $under+1; $i <= $over; $i++){
			$first +=  log($i);
		}
		# second sum
		for(my $i = 1; $i <= $over - $under; $i++){
			$second +=  log($i);
		}
		return $first-$second;
	}
}
