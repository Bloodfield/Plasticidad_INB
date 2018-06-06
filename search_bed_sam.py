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
#       $ python search_bed_sam.py [bed file] [sam file]
#       
#       outputs in stdout
#
#       Nota: el bed no debe tener cabezera
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#       MAIN
#--------------------------------------------------------------------

import sys

#constantes de programa
debugg = 1

#Ciclo de read para un archivo:
BedLine = []

if debugg : print("Generando array de busqueda <- BED")
with open(sys.argv[1]) as f:
    for line in f:
        if debugg : print("Line:",line)
        
        
        # Genera diccionario
        Splited = line.split()
        Search = Splited[3].split(";")[2]
        BedLine.append(Search)
        if debugg : print("BedLine append",Search)

if debugg : print("Buscando")
with open(sys.argv[2]) as f:
    for line in f:
        if (line[0]!="@"):
            Splited = line.split()
            if(BedLine.__contains__(Splited[0])):
                print(line)
