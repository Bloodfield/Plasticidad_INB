//	RMFilter in C

/*	TODO
 * 	Mejorar Caché
 * 	Hacer guía
 * 	Funcionalizal main
 */
#include <stdio.h>

#define Log_name "Log.txt"

long unsigned int overlap_range(long unsigned int Repeats_coordenates[][2],long unsigned int *flanco,int array_length);
FILE *Log;

int main(int argc, char *argv[]){
	
	Log = fopen(Log_name,"a");
	
	if (argc != 5){
		fprintf(Log,"Se tiene que escribir un archivo de referencia, con su número de líneas \n");
		fprintf(Log,"$ RMFilter Referencia.bin n_lineas n_bases porcentaje  \n");
		fprintf(Log,"$ porcentaje = porcentaje máximo de bases del flanco que corresponde a un repetido \n");
		fprintf(Log,"Ejemplo:\n");
		fprintf(Log,"\n");
		fprintf(Log,"$ RMFilter $HOME/rmsk_Crh1.bin 100 300 20 \n");
		fprintf(Log,"\n");
		fprintf(Log,"\n");
		fprintf(Log,"Este archivo de referencia se construye con el programa \"rmsk_files.sge\"\n");
		fprintf(Log,"Y la cantidad de lineas la encuentras en \"files_lines.txt\"\n");
		fprintf(Log,"\n");
		fprintf(Log,"\n");
		fclose(Log);
		return 1;
	}
	
	//	Obtiene argumentos
	
	char *RM_list_name;
	long unsigned int array_length = 0;
	long unsigned int qbases= 300;
	long unsigned int percentage= 20;
	fprintf(Log,"Name : %s \n",argv[1]);
	RM_list_name=argv[1];
	sscanf (argv[2],"%lu",&array_length);
	sscanf (argv[3],"%lu",&qbases);
	sscanf (argv[4],"%lu",&percentage);
	
	//	Obtiene la base de datos de Repeat Masker
	
	FILE *RM_table;
	RM_table = fopen(RM_list_name,"r");
	if(!RM_table){
		printf("Archivo : %s no existe\n",RM_list_name);
		return 1;
	}
	
	long unsigned int Repeats_coordenates[array_length][2];
	int index = 0;
	int x = 0, y=0;
	char c;
	
	while (fread(&Repeats_coordenates[x][y], sizeof(long unsigned int), 1, RM_table)== 1){
		index++;
		x= index/2;
		y=index%2;
	}
	fclose(RM_table);
	
	//	Revisa flancos
	long unsigned int eq_bases = qbases*percentage/100;
	char Chromosome[50];
	char Label[200];
	long unsigned int inicio, fin;
	while (!feof(stdin)){
		scanf("%s\t%lu\t%lu\t%[^\n]s\n",Chromosome, &inicio, &fin,Label);
		long unsigned int flanco[2];
		flanco[0]= inicio-qbases;
		flanco[1] = inicio;
		long unsigned int overlap = overlap_range(Repeats_coordenates,flanco,array_length);
		if(overlap < eq_bases){
			flanco[0]= fin;
			flanco[1] = fin + qbases;
			overlap = overlap_range(Repeats_coordenates,flanco,array_length);
			if(overlap < eq_bases){
				printf("%s\t%lu\t%lu\t%s\n",Chromosome, inicio, fin,Label);
			}
		}
	}
	
	fclose(Log);
	return 0;
}

//	Computa el overlap sobre "Repeats_coordenates" del intervalo en "Flanco"

long unsigned int overlap_range(long unsigned int Repeats_coordenates[][2],long unsigned int *flanco, int array_length){
	int prev_end = 0;
	int overlap = 0;
	int counter = 0;
	while(Repeats_coordenates[counter][1] < flanco[0] && counter < array_length){
		counter++;
	}
	prev_end = flanco[0];
	while(Repeats_coordenates[counter][0] < flanco[1] && counter < array_length){
		long unsigned int a = flanco[1] < Repeats_coordenates[counter][1] ? flanco[1] : Repeats_coordenates[counter][1];
		long unsigned int b = Repeats_coordenates[counter][0] > prev_end ? Repeats_coordenates[counter][0] : prev_end;
		overlap += a-b;
		prev_end = a;
		counter++;
	}
	return overlap;
}
