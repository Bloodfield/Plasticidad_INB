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
#       FUNCTIONS
#--------------------------------------------------------------------
#   1: Flag2Arr


#------------------------------
#       Flag2Arr
#------------------------------
#   Considera un numero entero que equivale a la bandera
#   de un SAM-file y regresa una lista que corresponde a
#   los campos que detectan las baderas del read

def Flag2Arr(Flag):
    FlagT=int(Flag)
    ret = [0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(12):
        ret[i] = FlagT % 2
        FlagT = int(FlagT / 2)
    return ret


#--------------------------------------------------------------------
#       MAIN
#--------------------------------------------------------------------

import re
import sys

#constantes de programa
GapLab = "DN"
MLab = "MSHP=X"
C_MIN = 20
G_MIN= 500

# CIGAR parce regexp
regex = re.compile("[0-9]+[M,N,D,S,H,P,=,X]")

with open(sys.argv[1]) as f:
    for line in f:
        if (line[0]!="@"):
            Splited = re.split('\W+',line)
            CIGAR = regex.findall(Splited[5])
            
                    
            print (CIGAR)
            print (Flag2Arr(Splited[1]))



