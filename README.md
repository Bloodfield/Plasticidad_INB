#	Plasticidad INB

##	

##	Descripción:

Estas herramientas estan hechas para poder llevar a cabo los pasos de detección y filtrado de deleciones somáticas en células de humano.
En lo particular toma los requerimientos en el laboratorio de neurodesarrollo en el [Instituto de Neurobiolgía (UNAM-Juriquilla)](http://www.inb.unam.mx/index.html)
De igual manera contiene scripts que autimatizan gran parte del proceso, desde la descarga de datos crudos, hasta su análisis de deleciones.

##	Consideraciones

La parte automática esta diseñada para ser ejecutada con el software Sun Grid Engine (SGE), con el cual se ingresan los trabajos de cómputo en el servidor del [Laboratorio de Visualización Científica y Avanzada (Lavis)](http://lavis.unam.mx/) de la UNAM-Juriquilla.

Es necesario revisar que se tienen a dispocición las siguientes dependencias:

+ SRA toolkit
+ fastx
+ fastqc
+ seqtk
+ Segemehl (0.11)
+ Samtools
+ [Bedtools](http://bedtools.readthedocs.io/)
+ fastc

##	Instalación

Primero hay que crear un folder $DIR/bin/ en tu carpeta de  en el cual se guardarán los programas.
La dirección DIR es la cual contendrá los datos as procesar por medio de los siguientes comandos (usando la dirección correspondiente).
La dirección contenida en $DIR será tomada como la capteta "HOME" en los programas *.sge que se crean con la intalación, por lo que tienen que elegirse con cuidado.
```
DIR="/mnt/Timina/mhernandez/echavez/"
mkdir $DIR/bin/
cd $DIR/bin/
```

Una vez que se tienen las dependencias instaladas y el archivo "$DIR/bin" creado, entoces se puede proceder a la descarga e instalación.
Para descargar primero hay que clonar el repositorio de github y entrar al folder respectivo:
(El clon puede ubicarse en  cualquier folder de tu usuario)
```
git clone https://github.com/Bloodfield/Plasticidad_INB.git
cd Plasticidad_INB
```

Dentro de la carpeta se encuentra el instalador "Install.sh"
Al ejecutarlo, se usará la direccón en "$DIR" por lo que se instalará en la capeta elegida.
Un solo clon del repositorio es capas de instalar los programas en las carpetas necesarias.
```
./Install $DIR
```
La ejecución instala los programas necesarios y copia los archivos *.sge automatizados

##	Referencias

+ La base de datos de repetidos se tomó del proyecto [RepeatMasker](http://www.repeatmasker.org/)
+ Las referencias de los centrómeros se obtubieron de la [anotación](https://www.ncbi.nlm.nih.gov/genome/guide/human/) del genoma humano de referencia hg38
+ el archivo de longitudes de cromosomas se obtubo de los recursos de [Bedtools](http://bedtools.readthedocs.io/)
##	TODO list

1. Mejorar proceso de instalación
2. Marcar referencias de bases de datos
3. Links de las dependencias
