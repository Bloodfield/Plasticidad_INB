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
int Del_Overlap(FILE *In_file,FILE *Log,unsigned overlap_th, unsigned score_th,char *Chr);
int add_line(FILE *In_file, unsigned *flanco_A, unsigned *flanco_B, unsigned *fin, FILE *Log);

int main(int argc, char *argv[]){
	
	FILE *Log;
	Log = fopen(Log_name,"a");
	
	char *In_file_name;
	unsigned score_th   = 1;
	unsigned overlap_th = 1;
	
	//	Entrada de datos
	
	if (argc != 5){
		err_message(Log);
		fclose(Log);
		return 1;
	}
	
	fprintf(Log,"Name : %s \n",argv[1]);
	In_file_name=argv[1];
	sscanf (argv[2],"%u",&score_th);
	sscanf (argv[3],"%u",&overlap_th);
	char *Chr = argv[4];
	
	if(overlap_th>100){
		fprintf(Log,"Porcentaje de overlap no válido \n");
		err_message(Log);
		fclose(Log);
		return 1;
	}
	
	// Ejecución del programa
	FILE *In_file;
	In_file= fopen(In_file_name,"r");
	if(!In_file){
		printf("Archivo : %s no existe\n",In_file_name);
		return 1;
	}
	unsigned Overflow = Del_Overlap(In_file,Log,overlap_th,score_th,Chr);
	
	if (Overflow==1){return 1;}
	
	// Cierre de archivos
	fclose(In_file);
	fclose(Log);
	
	return 0;
}

//	DeQueue el primer elemento y reordena la lista
unsigned de_queue(unsigned *array,unsigned end){
	unsigned temp=array[0];
	int i;
	for (i =0 ; i< end-1; i++)
		array[i]= array[i+1];
	array[end-1]=0;
	return temp;
}

//	Queue un elemento a la lista
unsigned queue(unsigned *array,unsigned element, unsigned end, FILE *Log){
	if(end + 1 >= Array_Size){
		fprintf(Log,"Err = Stack Overflow\n");
		return 1;
	}
	
	array[end] = element;
	
	return 0;
}

//	calcula el porcentaje de la referencia que es cubierta por el query
int percentage_overlap(unsigned reference , unsigned query ,unsigned *flanco_A,unsigned *flanco_B){
	unsigned FAR = flanco_A[reference];
	unsigned FBR = flanco_B[reference];
	unsigned FAQ = flanco_A[query];
	unsigned FBQ = flanco_B[query];
	if (FBQ < FAR || FAQ > FBR){return 0;}
	unsigned divisor = FBR - FAR;
	unsigned temp= 0;
	if (FAQ < FAR){
		if(FBQ < FBR){
			temp = (FBQ-FAR )*100;
		}else{
			return 100;
		}
	}else{
		if(FBQ < FBR){
			temp = (FBQ-FAQ )*100;
		}else{
			temp = (FBR-FAQ )*100;
		}
	}
	
	temp /= divisor;
	return temp;
}

int err_message(FILE *Log){
	fprintf(Log,"Se tiene que escribir un archivo bed para analizar \n");
	fprintf(Log,"El bed tiene que ser de un solo cromosoma, y estar ordenado de menor a mayor en base a la primera coordenada \n");
	fprintf(Log,"$ Del_Overlap [Name.bed] overlap_th score_th Chr\n");
	fprintf(Log,"Donde overlap_th es el porcentaje entero mínimo para decir que dos lecturas se sobrelapan y\n");
	fprintf(Log,"\t score_th es la mínima cantidad de sobrelapes para reportar una deleción\n");
	fprintf(Log,"\t Chr es el nombre del cromosoma que aparece en el bed de salida\n");
	fprintf(Log,"Emeplo:\n");
	fprintf(Log,"\n");
	fprintf(Log,"$ Del_Overlap $HOME/Deletions_Chr1.bed 50 5 Chr1\n");
	fprintf(Log,"\n");
	fprintf(Log,"El output se reporta en standard output \n");
	fprintf(Log,"\n");
	return 0;
}

int Del_Overlap(FILE *In_file,FILE *Log,unsigned overlap_th, unsigned score_th,char *Chr){
	
	unsigned flanco_A[Array_Size]={0};
	unsigned flanco_B[Array_Size]={0};
	unsigned score[Array_Size]={0};
	
	unsigned fin=0;
	
	//	Primer valor 
	if( add_line(In_file, flanco_A, flanco_B, &fin,Log)==1){
		fprintf(Log,"Something wrong with bed File: line 1\n");
		fclose(In_file);
		err_message(Log);
		fclose(Log);
	}
	
	//	Segundo valor
	if( add_line(In_file, flanco_A, flanco_B, &fin,Log)==1){
		fprintf(Log,"Something wrong with bed File: line 1\n");
		fclose(In_file);
		err_message(Log);
		fclose(Log);
	}
	
	if(flanco_A[fin-1] < flanco_A[fin-2]){
		fprintf(Log,"Bed file not sorted\n");
		err_message(Log);
		fclose(Log);
		return 1;
	}
	
	// 	printf("FA = %u FB = %u \n",FA,FB);
	
	int is_EOF = 0;
	unsigned Overflow = 0;
	
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
				}else if((Overflow = add_line(In_file, flanco_A, flanco_B, &fin,Log))==1){
					is_EOF = 1;
					FAX = FB0;
				}else{
					if(flanco_A[fin-1] < flanco_A[fin-2]){
						fprintf(Log,"Bed file not sorted\n");
						err_message(Log);
						fclose(Log);
						return 1;
					}
					FAX = flanco_A[fin-1];
				}
			}
		}
		
		FB0= flanco_B[0];
		FAX= flanco_A[1];
		
		unsigned FA0= flanco_A[0];
		unsigned FBX= flanco_B[1];
		
		unsigned FA = de_queue(flanco_A,fin);
		unsigned FB = de_queue(flanco_B,fin);
		unsigned temp_score = de_queue(score,fin);
		fin --;
		
		if (Overflow == 2){return 1;}
		
		if (temp_score >= score_th && FAX != FA0 && FBX != FB0){
			printf("%s\t%u\t%u\t%u\n",Chr,FA,FB,temp_score);
		}
		
		if (fin <= 1){
			Overflow = add_line(In_file, flanco_A, flanco_B, &fin,Log);
		}
		if (Overflow == 2){return 1;}
	}
	
	unsigned FA = de_queue(flanco_A,fin);
	unsigned FB = de_queue(flanco_B,fin);
	unsigned temp_score = de_queue(score,fin);
	fin --;
	
	if (temp_score > score_th){
		printf("%u\t%u\t%u\n",FA,FB,temp_score);
	}
	
	return 0;
	
}

int add_line(FILE *In_file, unsigned *flanco_A, unsigned *flanco_B, unsigned *fin, FILE *Log){
	
	char Dummy1[Str_len]={0};
	char Dummy2[Str_len]={0};
	unsigned FA,FB;
	unsigned Overflow= 0;
	
	//	Primer valor
	if(fscanf(In_file,"%s\t%u\t%u\t%[^\n]s\n",Dummy1,&FA,&FB,Dummy2)!= 4){
		
		return 1;
	}
	Overflow =  queue(flanco_A,FA,*fin,Log);
	Overflow += queue(flanco_B,FB,*fin,Log);
	(*fin)++;
	if (Overflow > 0){return  2;}
	
	return 0;
}
