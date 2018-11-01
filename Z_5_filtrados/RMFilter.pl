#!/usr/bin/perl

#--------------------------------------------------------------------
#       SAM_parce
#   Edgar Chávez Aparicio
#
#   El programa acepta un bedfile con las. deleciones a analizar
#	Revisa una cantidad definida de bases a ambos flancos
#	Posteriormente revisa que porcentaje de los flacos
#		empalma con RepeatMasker
#	Reporta solo los que cumplen con los parámetros definidos
#
#	__q_bases___
#	____________|DDDDDDDDDDDD|____________
#	flanco izq	GAP	  flanco der
#	
#	RRR_________		  ____RRRRRRRR		R=repetidos
#
#	MaxPercentage -> eq_bases de R permitidos en cada flanco
#
#	____________|DDDDDDDDDDDD|____________
#	flanco izq	GAP	  flanco der
#		:)			:(	-> no pasa
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#       Dependencias:
#       
#	rmsk.txt del organizmo a analizar, dividido en archivos por cromosoma
#		> Chr[Número o letra]*.txt
#	
#       Guía de uso:
#	Antes de usar, revisar que la configuración del programa es correcta
#		Línea 73
#       
#       $ perl SAM_parce_GAP.pl <--opts> "[param]" [bed file]
#
#	[bed file]	:	este archivo tiene que representar las coordenadas
#				de las deleciones, en un mismo cromosoma
#       
#       outputs in stdout
#
#	--Chromosome	Número del cromosoma de búsqueda
#			Default: 1
#	--max_per	Máximo porcentaje de sobrelape con RepeatMasker
#			Default: 20
#	--q_nbases	número de bases flanco que se usarán para la búsqueda
#			Default: 20
#
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#       TO_DO
#
#	Guía de uso
#	Prueba real
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#       MAIN
#--------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;
use List::Util qw[min max];

#	Parámetros de ejecución

# Defaults
my $Chromosome = 1;
my $MaxPercentage = 20.0;
my $q_bases = 300;


GetOptions	('Chromosome:i'=> \$Chromosome,
		'max_per:f' => \$MaxPercentage,
		'q_nbases:i' => \$q_bases);

#	Configuración del programa

#my $rmsk_dir = '/mnt/Timina/mhernandez/echavez/';
my $rmsk_dir = './';


#	Programa

my $eq_bases = $MaxPercentage*$q_bases/100;
$rmsk_dir= join "", $rmsk_dir, "rmsk_Chr", $Chromosome, ".txt";

while (my $line = <STDIN>) {
	
	## Flanco Izquierdo
	
	
# 	print("_Lines\n");
	
	open(my $fh_rmsk, '<:encoding(UTF-8)', $rmsk_dir)
	or die "Error con: '$rmsk_dir' $!";
	
	my $overlap = 0;
	my $coord = 0;
	my @repeat = split '\t', <$fh_rmsk>;
	my @deletion = split '\t', $line;
	my %flanco = (
		"init"  => $deletion[1]-$q_bases,
		"end" => $deletion[1]);
# 	print($flanco{init},"\t",$flanco{end},"\n");
	
	while($repeat[7] < $flanco{init} && ! eof($fh_rmsk)){
		@repeat = split '\t', <$fh_rmsk>;
# 		print($repeat[7],"\n");
	}
# 	print($repeat[6],"\t",$repeat[7],"\n");
	while($repeat[6] < $flanco{end} && ! eof($fh_rmsk)){
		if($repeat[7]>$flanco{end}){
			$repeat[7]=$flanco{end};
		}
		
		
		if ($repeat[6] < $flanco{init}){
			$overlap = $repeat[7] - $flanco{init};
			$coord = max($repeat[7],$coord);
		}elsif($repeat[6] > $coord){
			$overlap += $repeat[7] - $repeat[6];
			$coord = $repeat[7];
		}elsif($repeat[7] > $coord){
			$overlap += $repeat[7] - $coord;
			$coord = $repeat[7];
		}
		
		@repeat = split '\t', <$fh_rmsk>;
	}
	
	close($fh_rmsk);
	
	## Flanco Derecho
	
	if($overlap < $eq_bases){
	
		open(my $fh_rmsk, '<:encoding(UTF-8)', $rmsk_dir)
		or die "Error con: '$rmsk_dir' $!";
		
		$coord = 0;
		$overlap = 0;
		@repeat = split '\t', <$fh_rmsk>;
		$flanco{init}=$deletion[2];
		$flanco{end} =$deletion[2]+$q_bases;
		
# 		print($flanco{init},"\t",$flanco{end},"\n");
		
		while($repeat[7] < $flanco{init} && ! eof($fh_rmsk)){
			@repeat = split '\t', <$fh_rmsk>;
# 			print($repeat[7],"\n");
		}
# 		print($repeat[6],"\t",$repeat[7],"\n");
		while($repeat[6] < $flanco{end} && ! eof($fh_rmsk)){
			if($repeat[7]>$flanco{end}){
				$repeat[7]=$flanco{end};
			}
			
			
			if ($repeat[6] < $flanco{init}){
				$overlap = $repeat[7] - $flanco{init};
				$coord = max($repeat[7],$coord);
			}elsif($repeat[6] > $coord){
				$overlap += $repeat[7] - $repeat[6];
				$coord = $repeat[7];
			}elsif($repeat[7] > $coord){
				$overlap += $repeat[7] - $coord;
				$coord = $repeat[7];
			}
			
			@repeat = split '\t', <$fh_rmsk>;
		}
		
		close($fh_rmsk);
		
		##	Imprimir resultado
		
		if($overlap < $eq_bases){
			print($line);
		}
		
	}
	
	close($fh_rmsk);	
	
	
}

