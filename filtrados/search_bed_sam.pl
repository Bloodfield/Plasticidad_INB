#!/usr/bin/perl -n

#       search_bed_sam
#   Edgar Chávez Aparicio
#
#   Es un programa que permite buscar cada uno de los
#   Identificadores de un bed en un sam
#
#   Referencias:
#
#       Sequence Alignment/Map Format Specification
#           The SAM/BAM Format Specification Working Group
#           1 Jun 2017
#       Bedtools Manual
#           Aaron R. Quinlan and Ira M. Hall
#           University of Virginia
#           
#
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#       Guía de uso:
#       
#       $ perl search_bed_sam.py [bed file] [sam file]
#       
#       outputs in stdout
#
#       Nota: el bed no debe tener cabezera
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#       MAIN
#--------------------------------------------------------------------

use strict;
use warnings;
 
#	ABRIR ARCHIVO

# abre lista de BED
# my $file_BED = $ARGV[0];
# open(my $fh_BED, '<:encoding(UTF-8)', $file_BED)
#   or die "Error con: '$file_BED' $!"; 

my	@search = ();

while (my $line = <>) {
	# obtener datos
	chomp $line;
	# #   print "$line\n";
	my	@elems	= split '\t', $line;
	print "$elems[3]";
	my 	@split	= split ';' , $elems[3];
	@search=(@search,$split[2]);
}

close($fh_BED);

my	%test = map { $_, 1 } @search;

# Busca en el SAM
my $file_SAM = $ARGV[1];
open(my $fh_SAM, '<:encoding(UTF-8)', $file_SAM)
  or die "Error con: '$file_SAM' $!"; 

while (my $line = <$fh_SAM>) {
	if($line !~ m/^\@/){
		my	@elems	= split ' ', $line;
		if	( $test{ $elems[0] } ){
			print "$line"
		}
	}
}
close($fh_SAM);