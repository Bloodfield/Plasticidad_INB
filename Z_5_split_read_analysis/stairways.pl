#!/usr/bin/perl -n

##---------------------------------------------------------
#   DESCRIPCIÓN
#   
#   Este programa permite revisar los bloques de reads de un SAM
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
#   Se lee un parametro final de cad alinea del SAM
#       de la forma \[('[0-9]*[G,C]',)*'[0-9]*[G,C]'\]
#        Por ejemplo: ['58C', '1G', '42C']
#   
#   IN:	SAM ordenado
#
#       perl stairways.pl $ordenado.sam
#   
#   OUT (stdout):
#
#       Sam separado en bloques
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


#	FUNCIONES DE PROGRAMA

sub read_length{
	my @nums    =  ();
	my @sumand  =  ();
	my $sum     =  0;
	@nums       =  split /[^0-9]+/ ,$_[0];
	@sumand     =  grep /\S/, @nums;
# 	print "@sumand";
	$sum       +=  $_ for @sumand;
	return $sum
}
  
#	VARIABLES INTERNAS
#my @reads	= ();
my $count	= 0;	# $#reads
my $count_bloq	= 0;
my $last	= 0;
my $first	= 0;
my $f_min	= 0;
my $i_max	= 0;
my $i		= 0;
my $f		= 0;

#	LEER ARCHIVO

while (my $line = <>) {
	# obtener datos
	chomp $line;
	# #   print "$line\n";
	my @elems	= split '\t', $line;
	$i		= $elems[3];
	
	
	#   print "$elems[-1]\n";
	$f	= &read_length($elems[-1]);
	$f	+=$i-1;
# 	print "$f\t$i\n";
	#print "$row\n";
	
# 	print "$i\t$last\n";
	#	revisar bloque
	if ($i > $last){
		if($count_bloq > 0){
			my $l	=	$last - $first;
			my $dl	=	$i_max - $f_min;
			print "\@$count_bloq\t$count\t$l\t$dl\n";
			#print "\@$count_bloq\t$#reads\t$l\t$dl\n";
			#print join("\n",@reads),"\n";
			
		}
		#@reads	= ($line);
 		$count	= 1;
		$count_bloq++;
		$last	= $f;
		$first	= $i;
		$f_min	= $f;
		$i_max	= $i;
		
	}else{
		#@reads=(@reads,$line);
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
	print "$line\n";
  
}
my $l	=	$last - $first;
my $dl	=	$i_max - $f_min;
print "\@$count_bloq\t$count\t$l\t$dl\n";
#print "\@$count_bloq\t$#reads\t$l\t$dl\n";
#print join("\n",@reads),"\n";