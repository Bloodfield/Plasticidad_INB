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
#       Guía de uso:
#       
#       $ python SAM_parce.py [sam file]
#       
#       outputs in stdout
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
Clev_MIN = 1
Gap_MIN= 1
debugg = 0

# CIGAR parce regexp
regex = re.compile("[0-9]+[M,N,D,S,H,P,X,=]")

#Ciclo de read para un archivo:

with open(sys.argv[1]) as f:
    for line in f:
        
        
        #   Descarta comentarios
        if (line[0]!="@"):
            if debugg : print("LINE:",line)
            
            # separa campos
            Splited = line.split("\t")
            if debugg : print("Splited:",Splited)
            #obtiene CIGAR
            CIGAR = regex.findall(Splited[5])
            
            #Busqueda de los slplit reads en Gap_Form
            #   CIGAR -> Print (lineas que contengan split reads)
            
            if (len(CIGAR)>0):
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
                            Gap_Form.pop(size-i)
                        else:
                            Gap_Form[size-1-i]=Gap_Form[size-1-i][0:-1]+"G"
                    
                    #   Si es sun Cleverage
                    elif(Gap_Form[size-1-i][-1] in CIGAR_Tags["C"]):
                        if(Gap_Form[size-i][-1] == "C"):
                            Gap_Form[size-1-i]=str(int(Gap_Form[size-1-i][0:-1])+int(Gap_Form[size-i][0:-1]))+"C"
                            Gap_Form.pop(size-i)
                        else:
                            Gap_Form[size-1-i]=Gap_Form[size-1-i][0:-1]+"C"
                            
                    #   Si es un caracter invalido.
                    else:
                        print ("Invalid char in ", splited[0]," : ",Gap_Form[-1][-1])
                        
                ##  Parte 2 encontrar el patron _,"(>Clev_MIN)C","(>Gap_MIN)G","(>Clev_MIN)C",_
                
                S=0
                ### Maquina de estados para encontrar la secuencia
                for i in Gap_Form:
                    if debugg: print("_ S=",S,"i[0:-1]=",i[0:-1],"i[-1]",i[-1])
                    if (S==0):
                        if (int(i[0:-1])>=Clev_MIN and i[-1]=="C"):
                            S=1
                        else:
                            S=0
                    elif (S==1):
                        if (int(i[0:-1])>=Gap_MIN and i[-1]=="G"):
                            S=2
                        else:
                            S=0
                    elif (S==2):
                        if (int(i[0:-1])>=Clev_MIN and i[-1]=="C"):
                            S=3
                        else:
                            S=0
                    elif (S==3):
                        S=3
                    else:
                        S=0
                
                if(S==3):
                    print(line,"\t",Gap_Form)
                if debugg: print(Splited[0],"\t",CIGAR,"\t",Gap_Form)
            else: print(Splited[0],"Not Valid")
print("Done")