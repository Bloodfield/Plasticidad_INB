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
CIGAR_Tags = {"G":"DN","C":"MSHP=X"}
Clev_MIN = 20
Gap_MIN= 500
debugg = 0

# CIGAR parce regexp
regex = re.compile("[0-9]+[M,N,D,S,H,P,=,X]")

#Ciclo de read para un archivo:

with open(sys.argv[1]) as f:
    for line in f:
        
        #   Descarta comentarios
        if (line[0]!="@"):
            
            # separa campos
            Splited = re.split('\W+',line)
            
            #obtiene CIGAR
            CIGAR = regex.findall(Splited[5])
            
            #Busqueda de los slplit reads en Gap_Form
            #   CIGAR -> Print (lineas que contengan split reads)
            
            ##  Parte 1 encontrar los Gaps y los Clevareges
            ##      CIGAR -> Gap_Form
            Gap_Form=CIGAR.copy()
            if(Gap_Form[-1][-1] in CIGAR_Tags["G"]):
                Gap_Form[-1]=Gap_Form[-1][0:-1]+"G"
            elif(Gap_Form[-1][-1] in CIGAR_Tags["C"]):
                Gap_Form[-1]=Gap_Form[-1][0:-1]+"C"
            else:
                print ("Invalid char in ", splited[0]," : ",Gap_Form[-1][-1])
            if debugg : print("_",Gap_Form)
            size = len(Gap_Form)-1
            for i in range(size):
                if debugg : print("_ _",i," ",Gap_Form)
                #   Si es sun gap
                if(Gap_Form[size-1-i][-1] in CIGAR_Tags["G"]):
                    if(Gap_Form[size-i][-1] == "G"):
                        Gap_Form[size-1-i]=str(int(Gap_Form[size-1-i][0:-1])+int(Gap_Form[size-i][0:-1]))+"G"
                        Gap_Form.remove(size-i)
                    else:
                        Gap_Form[size-1-i]=Gap_Form[size-1-i][0:-1]+"G"
                
                #   Si es sun Cleverage
                elif(Gap_Form[size-1-i][-1] in CIGAR_Tags["C"]):
                    if(Gap_Form[size-i][-1] == "C"):
                        Gap_Form[size-1-i]=str(int(Gap_Form[size-1-i][0:-1])+int(Gap_Form[size-i][0:-1]))+"C"
                        Gap_Form.remove(size-i)
                    else:
                        Gap_Form[size-1-i]=Gap_Form[size-1-i][0:-1]+"C"
                        
                #   Si es un caracter invalido.
                else:
                    print ("Invalid char in ", splited[0]," : ",Gap_Form[-1][-1])
                    
            ##  Parte 2 encontrar el patron _,"(>Clev_MIN)C","(>Gap_MIN)G","(>Clev_MIN)C",_
            
            
            print(Splited[0],"\t",CIGAR,"\t",Gap_Form)
                    
                    
            


