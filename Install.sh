#!/bin/bash

#	Revisar la entrada de argumentos
if [ $# != 2 ]
then
	echo "Error : Arguments where not properly assigned"
	echo "	or DIR has forbidden caracters (blank spaces for example)"
	exit 1
fi

#	Lectura del folder de referencia

HOME="$1"
email="$2"

if [ ! -d $HOME/bin ]
then
	echo "Error : Install file not found"
	echo "	$HOME/bin"
	exit 1
fi

echo "Instalando / Actualizando en $HOME"

echo "	Instalar los Programas de C"

gcc -o $HOME/bin/SAM_parce_GAP Z_5_filtrados/SAM_parce_GAP.c
gcc -o $HOME/bin/RMFilter Z_5_filtrados/RMFilter.c
gcc -o $HOME/bin/Del_Overlap Z_5_filtrados/Del_Overlap.c
gcc -o $HOME/bin/Allelic_Filter Z_6_Comparaciones/Allelic_Filter.c
gcc -o $HOME/bin/Coverage_count Z_5_filtrados/Coverage_count.c

echo "	Obtener base de datos de Repeat Masker"

if [ ! -d $HOME/bin/rmsk ]
then
	mkdir $HOME/bin/rmsk
fi
tar -xzf Z_5_filtrados/rmsk.tar.gz -C $HOME/bin/rmsk/

echo "	Copiar otras utilidades"

cp -f Z_6_Comparaciones/bed2hist.r $HOME/bin/bed2hist.r

cp -f Z_6_Comparaciones/human.hg38.genome $HOME/bin/human.hg38.genome
cp -f Z_6_Comparaciones/centromere_coords.bed $HOME/bin/centromere_coords.bed
cp -f Z_6_Comparaciones/genes_gdc_small_ucsc.bed $HOME/bin/genes_gdc_small_ucsc.bed

echo "	Copiar sge's en la carpeta HOME"

cp Z_1_DataBase_Genomes/SRA/sra_par_download.sge $HOME/sra_pipe_download.sge
cp Z_5_filtrados/Circular.sge $HOME/Filter.sge
cp Z_6_Comparaciones/BED_histogram.sge $HOME/BED_histogram.sge
cp Z_6_Comparaciones/SAM_view_set.sge $HOME/SAM_view_set.sge  
cp Z_6_Comparaciones/Recurrentes.sge $HOME/Recurrentes.sge

DIR=$(echo $HOME | sed 's/\//\\\//g')
sed -i "s/^HOME=\"\/mnt\/Timina\/mhernandez\/echavez\"$/HOME=\"$DIR\"/g" $HOME/sra_pipe_download.sge
sed -i "s/^HOME=\"\/mnt\/Timina\/mhernandez\/echavez\"$/HOME=\"$DIR\"/g" $HOME/Filter.sge
sed -i "s/^HOME=\"\/mnt\/Timina\/mhernandez\/echavez\"$/HOME=\"$DIR\"/g" $HOME/BED_histogram.sge
sed -i "s/^HOME=\"\/mnt\/Timina\/mhernandez\/echavez\"$/HOME=\"$DIR\"/g" $HOME/SAM_view_set.sge  
sed -i "s/^HOME=\"\/mnt\/Timina\/mhernandez\/echavez\"$/HOME=\"$DIR\"/g" $HOME/Recurrentes.sge

sed -i "s/^#\$ -M echavezaparicio@gmail.com/#\$ -M ${email}/g" $HOME/sra_pipe_download.sge
sed -i "s/^email=\"echavezaparicio@gmail.com\"/email=\"${email}\"/g" $HOME/sra_pipe_download.sge
sed -i "s/^#\$ -M echavezaparicio@gmail.com/#\$ -M ${email}/g" $HOME/Filter.sge
sed -i "s/^email=\"echavezaparicio@gmail.com\"/email=\"${email}\"/g" $HOME/Filter.sge
sed -i "s/^#\$ -M echavezaparicio@gmail.com/#\$ -M ${email}/g" $HOME/BED_histogram.sge
sed -i "s/^email=\"echavezaparicio@gmail.com\"/email=\"${email}\"/g" $HOME/BED_histogram.sge
sed -i "s/^#\$ -M echavezaparicio@gmail.com/#\$ -M ${email}/g" $HOME/SAM_view_set.sge  
sed -i "s/^email=\"echavezaparicio@gmail.com\"/email=\"${email}\"/g" $HOME/SAM_view_set.sge  
sed -i "s/^#\$ -M echavezaparicio@gmail.com/#\$ -M ${email}/g" $HOME/Recurrentes.sge
sed -i "s/^email=\"echavezaparicio@gmail.com\"/email=\"${email}\"/g" $HOME/Recurrentes.sge

echo "Terminado"
