#!/usr/bin/perl

##---------------------------------------------------------
#   DESCRIPCIÓN
#   
#   Este programa permite revisar los bloques de reads de un BED
#   Un bloque son reads que sobrelapan
#
#   ________________
#       ________________
#     ________________
#           _______
#       ___         _____
#                              ________________
#                                    ________________
#
#        Bloque 1                   Bloque 2
#
#   
#   
#   Importante:
#   
#   IN:	BED ordenado
#
#       perl stairways_BED.pl [$ordenado.bed]
#   
#   OUT (stdout):
#
#       BED separado en bloques
#       
##---------------------------------------------------------

##---------------------------------------------------------
#   TO DO LIST
#   
#
#   Guía de uso
##---------------------------------------------------------


use strict;
use warnings;
 
#	ABRIR ARCHIVO

# my $filename = $ARGV[0];
# open(my $fh, '<:encoding(UTF-8)', $filename)
#   or die "Error con: '$filename' $!"; 
  
#	VARIABLES INTERNAS

my $count	= 0;	# $#reads
my $count_bloq	= 0;
my $last	= 0;
my $first	= 0;
my $f_min	= 0;
my $i_max	= 0;
my $i		= 0;
my $f		= 0;
my $l_min	= 0;
my $l_avr	= 0.0;
my $l_max	= 0;

#	LEER ARCHIVO

while (my $line = <STDIN>) {
	
	# obtener datos
	chomp $line;
	my @elems	= split '\t', $line;
	$i		= $elems[1];
	$f		= $elems[2];
	
	#	revisar bloque
	if ($i > $last){
		
		#	Imprimir separador
		if($count_bloq > 0){
			my $l	=	$last - $first;
			my $dl	=	$i_max - $f_min;
			print "\@${first} _____\t${i_max}\t${f_min}\t____ ${last}\n";
			print "\@$count_bloq,n=$count,l=$l,d=$dl\n";
			print "\@Distr=[$l_min,$l_avr,$l_max]\n";
		}
		
		#	RESET
		$count	= 1;
		$count_bloq++;
		$last	= $f;
		$first	= $i;
		$f_min	= $f;
		$i_max	= $i;
		$l_max	= $f-$i;
		$l_min	= $f-$i;
		$l_avr	= 0;
		
	}
	
	#	Añadir read a bloque
	else{
 		$count++;
		if($f> $last){
			$last=$f;
		}
		if($f_min > $f){
			$f_min=$f;
		}
		if($i_max < $i){
			$i_max=$i;
		}
	}
	
	#	Length Stats
	my $l = $f-$i;
	if($l>$l_max){$l_max=$l;}
	elsif($l<$l_min){$l_min=$l;}
	$l_avr = (($count-1)*$l_avr+$l)/$count;
	
	print "$line\n";
  
}
my $l	=	$last - $first;
my $dl	=	$i_max - $f_min;
print "\@${first} _____\t${i_max}\t${f_min}\t____ ${last}\n";
print "\@$count_bloq,n=$count,l=$l,d=$dl\n";
print "\@Distr=[$l_min,$l_avr,$l_max]\n";