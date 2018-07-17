#!/usr/bin/perl

#--------------------------------------------------------------------
#       SAM_parce
#   Edgar Chávez Aparicio
#
#   Es un programa que permite la lectura y análisis de SAM-files
#	Encuetra Gaps con los parametros indicados.
#	Se reportan en formato BED
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
#       $ perl SAM_parce_GAP.pl [sam file]
#       
#       outputs in stdout
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#       MAIN
#--------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;

my	$C_min	= 1;
my	$G_min	= 1;
my	$G_max	= 100000;

#	ABRIR ARCHIVO

# abre SAM
# my $file_SAM = $ARGV[0];
# open(my $fh_SAM, '<:encoding(UTF-8)', $file_SAM)
#   or die "Error con: '$file_SAM' $!"; 

#	Prin Head
#print "track name=\"GAP_List\" description=\"Gaps of $G_min bases between $C_min alignments\" \n";

GetOptions	('MinClev:i'=> \$C_min,
		'MaxGap:i'=> \$G_max,
		'MinGap:i'=> \$G_min);


while (my $line = <STDIN>) {
	
	if($line !~ m/^\@/){
		
		#	Estract info and translate CIGAR
		my	@elems	= split '\t', $line;
		$elems[5] =~ tr/[DN]/G/;	#GAP
		$elems[5] =~ tr/[SHP]/O/;	#null
		$elems[5] =~ tr/[M=X]/C/;	#Cleverage
		my	@CIGAR = ( $elems[5] =~ m/([0-9]+[GC])/g );
		
		
		#	Máquina de estados
		
		my	$state	= "S0";
		my	$C1	= 0;
		my	$DC2	= 0;
		my	$C2	= 0;
		my	$G	= 0;
		my	$inicio	= 0;
		my	$fin	= 0;
		
		for my $i (0 .. $#CIGAR){
			if ($CIGAR[$i] =~ m/C$/){
				$C2 +=  substr $CIGAR[$i], 0, -1;
			}
		}
		
		for my $i (0 .. $#CIGAR){
		
			#	Stand By
			if ($state eq "S0") {
				if($CIGAR[$i] =~ m/C$/){
					$C1 +=  substr $CIGAR[$i], 0, -1;
					$state = "C1";
				}
			}
			
			#	Clevarage 1
			elsif ($state eq "C1") {
				if($CIGAR[$i] =~ m/C$/){
					$C1 +=  substr $CIGAR[$i], 0, -1;
				}elsif($CIGAR[$i] =~ m/G$/){
					$G +=   substr $CIGAR[$i], 0, -1;
					$state = "G";
				}
			}
			
			#	Gap
			elsif ($state eq "G") {
				if($CIGAR[$i] =~ m/G$/){
					$G +=  substr $CIGAR[$i], 0, -1;
				}elsif($CIGAR[$i] =~ m/C$/){
					$DC2 +=   substr $CIGAR[$i], 0, -1;
					$state = "C2";
					$C2 -= $C1;
				}
			}
			
			#	Clevarage 2
			elsif ($state eq "C2") {
				if($CIGAR[$i] =~ m/C$/){
					$DC2 +=  substr $CIGAR[$i -1], 0, -1;
				}elsif($CIGAR[$i] =~ m/G$/){
				
					#	Test -> print
					if($C1 >= $C_min && $C2 >= $C_min && $G >= $G_min && $G <= $G_max){
						$inicio	= $elems[3] + $C1;
						$fin	= $inicio + $G;
						print "$elems[2]\t$inicio\t$fin\t$elems[0]\n";
					}
					
					#	RESET
					$C1	+= $DC2;
					$DC2	= 0;
					$G	= substr $CIGAR[$i], 0, -1;
					$inicio	= 0;
					$fin	= 0;
					$state	= "C1";
					
				}
			}
			
			#	FSM Err
			else {print "Err if case";}
		}
		
		#	Last Step
		if($C1 >= $C_min && $C2 >= $C_min && $G >= $G_min && $G <= $G_max){
			$inicio	= $elems[3] + $C1;
			$fin	= $inicio + $G;
			print "$elems[2]\t$inicio\t$fin\t$elems[0]\n";
		}
	}
}