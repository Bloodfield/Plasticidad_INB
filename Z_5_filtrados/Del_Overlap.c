//	RMFilter in C

/*	Descripción
 * 	Lee un Bed file ordenado por la primera coordenada
 * 	Imprime las lecturas que estén sobre un umbral de:
 * 		Sobrelape
 * 		De frecuencia en formato "coordenada_inicio \t coordenada_final \t cuenta de frecuencia
 * 
 */

/*	TODO
 * 	Mejorar Caché
 * 	Hacer guía
 * 
 */
#include <stdio.h>
#include <string.h>

#define Log_name "Log.txt"
#define Array_Size 300
#define Str_len 300

unsigned de_queue(unsigned *array,unsigned end);
unsigned queue(unsigned *array,unsigned element, unsigned end, FILE *Log);
int percentage_overlap(unsigned reference , unsigned query ,unsigned *flanco_A,unsigned *flanco_B);
int err_message(FILE *Log);

int main(int argc, char *argv[]){
	
	FILE *Log;
	Log = fopen(Log_name,"a");
	
	char *In_file_name;
	unsigned score_th   = 1;
	unsigned overlap_th = 1;
	
	// 	fprintf(Log,"Key : %s \n",argv[1]);
	
	if (argc != 4){
		err_message(Log);
		fclose(Log);
		return 1;
	}
	
	fprintf(Log,"Name : %s \n",argv[1]);
	In_file_name=argv[1];
	sscanf (argv[2],"%u",&score_th);
	sscanf (argv[3],"%u",&overlap_th);
	
	if(overlap_th>100){
		fprintf(Log,"Porcentaje de overlap no válido \n");
		err_message(Log);
		fclose(Log);
		return 1;
	}
	
	unsigned flanco_A[Array_Size]={0};
	unsigned flanco_B[Array_Size]={0};
	unsigned score[Array_Size]={0};
	
	unsigned fin=0;
	
	FILE *In_file;
	In_file= fopen(In_file_name,"r");
	char Dummy1[Str_len]={0};
	char Dummy2[Str_len]={0};
	unsigned FA,FB;
	
	//	Primer valor
	if(fscanf(In_file,"%s\t%u\t%u\t%[^\n]s\n",Dummy1,&FA,&FB,Dummy2)!= 4){
		fprintf(Log,"Something wrong with bed File: line 1\n");
// 		int test =fscanf(In_file,"%s\t%u\t%u\t%[^\n]s\n",Dummy1,&FA,&FB,Dummy2);
// 		fprintf(Log,"Next: %d\n",test);
		fprintf(Log,"Dummy1 = %s FA = %u FB = %u Dummy2 = %s\n",Dummy1,FA,FB,Dummy2);
		fclose(In_file);
		err_message(Log);
		fclose(Log);
		return 1;
	}
	queue(flanco_A,FA,fin,Log);
	queue(flanco_B,FB,fin,Log);
	fin++;
	
// 	printf("FA = %u FB = %u \n",FA,FB);
	
	//	Segundo valor
	if(fscanf(In_file,"%s\t%u\t%u\t%[^\n]s\n",Dummy1,&FA,&FB,Dummy2)!= 4){
		fprintf(Log,"Something wrong with bed File: line 2\n");
		fprintf(Log,"Dummy1 = %s FA = %u FB = %u Dummy2 = %s\n",Dummy1,FA,FB,Dummy2);
		fclose(In_file);
		err_message(Log);
		fclose(Log);
		return 1;
	}
	queue(flanco_A,FA,fin,Log);
	queue(flanco_B,FB,fin,Log);
	fin++;
	
// 	printf("FA = %u FB = %u \n",FA,FB);
	
	int is_EOF = 0;
	
	while(fin > 1){
		int X = 1;
		unsigned FB0= flanco_B[0];
		unsigned FAX= flanco_A[X];
		score[0]++;
		while(FB0>FAX){
			int perA0 = percentage_overlap(X,0,flanco_A,flanco_B);
			int per0A = percentage_overlap(0,X,flanco_A,flanco_B);
			
			if(per0A > overlap_th && perA0 > overlap_th){
				score[0]++;
				score[X]++;
			}
			X++;
			if(X>=fin){
				if(is_EOF){
					FAX = FB0;
				}else if(fscanf(In_file,"%s\t%u\t%u\t%[^\n]s\n",Dummy1,&FA,&FB,Dummy2)!= 4){
					is_EOF = 1;
					FAX = FB0;
				}else{
					queue(flanco_A,FA,fin,Log);
					queue(flanco_B,FB,fin,Log);
					fin++;
					FAX = FA;
				}
			}
		}
		FA = de_queue(flanco_A,fin);
		FB = de_queue(flanco_B,fin);
		unsigned temp_score = de_queue(score,fin);
		fin --;
		
		if (temp_score >= score_th){
			printf("%u\t%u\t%u\n",FA,FB,temp_score);
		}
	}
	
	FA = de_queue(flanco_A,fin);
	FB = de_queue(flanco_B,fin);
	unsigned temp_score = de_queue(score,fin);
	fin --;
	
	if (temp_score > score_th){
		printf("%u\t%u\t%u\n",FA,FB,temp_score);
	}
	
	fclose(In_file);
	fclose(Log);
	
	return 0;
}

unsigned de_queue(unsigned *array,unsigned end){
	unsigned temp=array[0];
	int i;
	for (i =0 ; i< end-1; i++)
		array[i]= array[i+1];
	array[end-1]=0;
	return temp;
}

unsigned queue(unsigned *array,unsigned element, unsigned end, FILE *Log){
	if(end + 1 >= Array_Size){
		fprintf(Log,"Stack Overflow");
		return 1;
	}
	
	array[end] = element;
	
	return 0;
}

int percentage_overlap(unsigned reference , unsigned query ,unsigned *flanco_A,unsigned *flanco_B){
	
	return 100;
}

int err_message(FILE *Log){
	fprintf(Log,"Se tiene que escribir un archivo bed para analizar \n");
	fprintf(Log,"El bed tiene que ser de un solo cromosoma, y estar ordenado de menor a mayor en base a la primera coordenada \n");
	fprintf(Log,"$ Del_Overlap [Name.bed] overlap_th score_th \n");
	fprintf(Log,"Donde overlap_th es el porcentaje entero mínimo para decir que dos lecturas se sobrelapan y\n");
	fprintf(Log,"\t score_th es la mínima cantidad de sobrelapes para reportar una deleción\n");
	fprintf(Log,"Emeplo:\n");
	fprintf(Log,"\n");
	fprintf(Log,"$ Del_Overlap $HOME/Deletions_Crh1.bed 50 5 \n");
	fprintf(Log,"\n");
	fprintf(Log,"El output se reporta en standard output \n");
	fprintf(Log,"\n");
	return 0;
}
