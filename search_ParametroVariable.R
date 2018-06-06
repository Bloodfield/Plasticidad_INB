#   Script para la búsqueda de parámetros variables de NCBI
#   Por: Chávez Aparicio Edgar Iván

#   Fecha: 16-Abril-2018

#   Dependencias:
#       rentrez

#   INPUT

BasicQuery <- "Human[ORGN] AND strategy wgs[FILT] AND Public [ACS] AND genomic [FILT]"
Base_Datos <- "sra"

ParametroVariable <- "[MBS]"
#Rango <- 360000:500000
#Rango <- append(225000:369999,500000:550000)
Rango <- 800:1200

UID_size <- 10
Repeats <- 50

#   OUTPUT
#   Lista de tamaño por valor del parámetro (No despliega resultados 0)
Exit_report <- "Report_Param.txt"
#   Id's de resultado al retmax especificado
Exit_UIDs <- "UIDs.txt"

#---------

#   Librería Entrez utilities
library("rentrez") 

# Script

# creación de arreglos

Longitud <- length(Rango)
Dividendo <- length(Rango)/Repeats
Residuo <- length(Rango)%%Repeats
Search_vals <- array(dim=c(Dividendo+1,Repeats))

for(i in 1:(Dividendo)){
    Search_vals[i,] <- Rango[(((i-1)*Repeats)+1):((i)*Repeats)]
}
Search_vals[Dividendo+1,c(1:(Residuo))] <- Rango[(Longitud-Residuo+1):Longitud]

#   Búsqueda

report_file <- file(Exit_report,"w")
UIDs_file <- file(Exit_UIDs,"w")

#Cabeceras de los archivos de salida
text <- paste(ParametroVariable,"\tCount")
write(text,report_file,ncolumns= 2)

write("UIDs",UIDs_file,ncolumns= 1)

print("Init")

for (i in 1:length(Search_vals[,1])){
    
    valor <-Search_vals[i,]
    
    #   Busca en la base de datos
    
    query <- paste( "OR",valor[2:length(valor)],rep(ParametroVariable,Repeats-1),sep=" ",collapse=" ")
    query <- paste(valor[1],ParametroVariable,query)
    query <- paste("(",BasicQuery, ") AND (", query, ")")
    
    search_res <- entrez_search(db=Base_Datos, term=query, retmax=UID_size)
    
    print(paste(query,search_res$count))
    
    if(search_res$count>0){
        #   Guarda los datos Exit_report
        text <- paste(valor[1],search_res$count,sep="\t")
        write(text,report_file,ncolumns= 1)
        
        #   Guarda los datos Exit_UIDs
        write(search_res$ids,UIDs_file,ncolumns= 1)
    }
    
}

close(report_file)
close(UIDs_file)

print("Exit")

# To Do List