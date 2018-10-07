#!/usr/bin/perl -n

#--------------------------------------------------------------------
#       SAM_parce
#   Edgar Chávez Aparicio
#
#   Es un programa que permite la lectura y análisis de SAM-files
#
#   Referencias:
#
#       Sequence Alignment/Map Format Specification
#           The SAM/BAM Format Specification Working Group
#           1 Jun 2017
#
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#       Guía de uso:
#       
#       $ perl SAM_parce.pl [sam file]
#       
#       outputs in stdout
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#       MAIN
#--------------------------------------------------------------------

use strict;
use warnings;
 

my	$C_min	= 1;
my	$G_min	= 1;

#	ABRIR ARCHIVO

# abre SAM
# my $file_SAM = $ARGV[0];
# open(my $fh_SAM, '<:encoding(UTF-8)', $file_SAM)
#   or die "Error con: '$file_SAM' $!"; 

while (my $line = <>) {
	if($line !~ m/^\@/){
# 		print "$line";
		my	@elems	= split '\t', $line;
		$elems[5] =~ tr/[DN]/G/;
		$elems[5] =~ tr/[MSHP=X]/C/;
		my	@CIGAR = ( $elems[5] =~ m/([0-9]+[GC])/g );	#not working
		
		
		#	Máquina de estados
		
		my	$state	= "S0";
		my	$C1	= 0;
		my	$C2	= 0;
		my	$G	= 0;
		my	$flag	= 0;
		
		for my $i (0 .. $#CIGAR){
			if ($state eq "S0") {
				if($CIGAR[$i] =~ m/C$/){
					$C1 +=  substr $CIGAR[$i], 0, -1;
					$state = "S1";
				}
			}
			elsif ($state eq "S1") {
				if($CIGAR[$i] =~ m/C$/){
					$C1 +=  substr $CIGAR[$i], 0, -1;
				}elsif($CIGAR[$i] =~ m/G$/){
					$G +=   substr $CIGAR[$i], 0, -1;
					$state = "S2";
				}else{print "Err in CIGAR\n";}
			}
			elsif ($state eq "S2") {
				if($CIGAR[$i] =~ m/G$/){
					$G +=  substr $CIGAR[$i], 0, -1;
				}elsif($CIGAR[$i] =~ m/C$/){
					$C2 +=   substr $CIGAR[$i], 0, -1;
					$state = "S3";
				}else{print "Err in CIGAR\n";}
			}
			elsif ($state eq "S3") {
				if($CIGAR[$i] =~ m/C$/){
					$C1 +=  substr $CIGAR[$i -1], 0, -1;
				}elsif($CIGAR[$i] =~ m/G$/){
					if($C1 > $C_min && $C2 > $C_min && $G > $G_min){
						$flag = 1;
					}
					$C1 = $C2;
					$C2 = 0;
					$G =   substr $CIGAR[$i], 0, -1;
					$state = "S2";
				}else{print "Err in CIGAR\n";}
			}
			else {print "Err if case";}
		}
		if($C1 > $C_min && $C2 > $C_min && $G > $G_min){$flag = 1;}
		
		
		if($flag == 1) {print "$line";}
		
		
	}
}