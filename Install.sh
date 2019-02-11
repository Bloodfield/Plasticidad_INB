#!/bin/bash

#	Configuración

HOME="/mnt/Timina/mhernandez/echavez/"


echo "Installing / Updating ..."

#	Instalar los Programas de C

gcc -o ../SAM_parce_GAP Z_5_filtrados/SAM_parce_GAP.c
gcc -o ../RMFilter Z_5_filtrados/RMFilter.c
gcc -o ../Del_Overlap Z_5_filtrados/Del_Overlap.c
gcc -o ../Allelic_Filter Z_6_Comparaciones/Allelic_Filter.c
gcc -o ../Coverage_count Z_5_filtrados/Coverage_count.c
#gcc -o ../RM_txt2bin Z_5_filtrados/RM_txt2bin.c

#	Obtener base de datos de Repeat Masker

#	

#	copiar sge's en la carpeta HOME

cp Z_1_DataBase_Genomes/SRA/sra_par_download.sge $HOME/sra_pipe_download.sge
cp Z_5_filtrados/Circular.sge $HOME/Filter.sge
cp Z_6_Comparaciones/BED_histogram.sge $HOME/BED_histogram.sge
cp Z_6_Comparaciones/bed2hist.r $HOME/bin/bed2hist.r
cp Z_6_Comparaciones/centromere_coords.bed $HOME/bin/centromere_coords.bed
co Z_6_Comparaciones/SAM_view_set.sge $HOME/SAM_view_set.sge  

echo Finish
