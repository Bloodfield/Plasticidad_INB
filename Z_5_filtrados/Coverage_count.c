#include <stdio.h>

#define MAX_ARR 500000
#define Log_name "Log.txt"
#define Str_len	40
#define n_chr	460

//	Sistema
int print_help();

//	Algoritmo
int Coverage_Count();
int fill_buffer(char *Chr_read,char *Chr_last, int inicio,int *Line_coord_1, int *Line_coord_2,int *size);
int count_score(int *Line_coord_1, int *Line_coord_2,int *Count,int limit,int inicio,int fin,int cover_score);
int print_buffer(int *Line_coord_1, int *Line_coord_2,int *Count,int limit,int *size,char *Chr_read,int fin);

//	Strings
int str_copy(char *Str1,char *Str2);
int str_clear(char *Str);
int str_eq(char *Str1,char *Str2);


//	Basic numerics
int min(int a, int b);
int max(int a, int b);

//	Arrays

int recorrer_array(int *array , int size, int n);

//	Globales
FILE *Log;
FILE *bed2bam_fh;

int main(int argc, char *argv[]){
	Log = fopen(Log_name,"a");
	//	Input 
	char *bed2bam;
	
	if (argc != 2){
		fprintf(Log,"argc = %d\n",argc);
		print_help();
		fclose(Log);
		return 1;
	}
	
	fprintf(Log,"Name_bed2bam : %s \n",argv[1]);
	bed2bam=argv[1];
	
	//	Abrir y probar archivos
	
	
	bed2bam_fh= fopen(bed2bam,"r");
	
	int Watchdog =0;
	if(!bed2bam_fh ){
		fprintf(Log,"Archivo : %s no existe\n",bed2bam);
		fclose(bed2bam_fh);
		fclose(Log);
		return 1;
	}
	
	//	Programa
	Watchdog = Coverage_Count();
	
	if(Watchdog){
		fprintf(Log,"Failled algorithm\n");
		fclose(bed2bam_fh);
		fclose(Log);
		return 1;
	}
	
	fclose(bed2bam_fh);
	fclose(Log);
	return 0;
	
}

int print_help(){
	
	fprintf(Log,"Error. Mostrando ayuda: \n");
	fprintf(Log," \n");
	fprintf(Log,"Coverage_Count permite realizar la cuenta de las coordenadas leidas en BED2BAM \n");
	fprintf(Log," Utiliza la saida de \"bedtools genomecov\" en formato BED graph\n");
	fprintf(Log,"Se tiene que tener el mismo orden de cromosomas en ambos. \n");
	fprintf(Log," Las coordenadas tienen que estar ordenadas igual\n");
	fprintf(Log," \n");
	fprintf(Log,"Su uso es: \n");
	fprintf(Log," \n");
	fprintf(Log,"$ Coverage_Count Archivo_BED2Bam.bed \n");
	fprintf(Log," Recibiendo por stdin\n");
	fprintf(Log," \n");
	fprintf(Log," \n");
	fprintf(Log,"Ejemplo:\n");
	fprintf(Log,"\n");
	fprintf(Log,"$ bedtools genomecov -split -bg -ibam SRR822853_sorted.bam  \\ \n");
	fprintf(Log,"\t | overage_Count SRR822853_BED2Bam.bed \\ \n");
	fprintf(Log,"\t > SRR822853_bedcov.bed \n");
	fprintf(Log,"\n");
	return 0;
}

int Coverage_Count(){
	int Line_coord_1[MAX_ARR]={0};
	int Line_coord_2[MAX_ARR]={0};
	int Count[MAX_ARR]={0};
	int Watchdog = 0;
	int size = 0;
	char Chr_last[Str_len]={0};
	
	while( !feof(stdin)){
		
		char Chr_read[Str_len]={0};
		int inicio = 0;
		int fin = 0;
		int cover_score = 0 ;
		
		//	Get line
		if(scanf("%s\t%u\t%u\t%u\n",Chr_read,&inicio,&fin,&cover_score)!= 4){
			fprintf(Log,"Archivo : Not a BED graph format in the stdin file\n");
			return 1;
		}
		
		Watchdog = fill_buffer(Chr_read,Chr_last,fin,Line_coord_1, Line_coord_2,&size);
		if(Watchdog){
			fprintf(Log,"Failled fill buffer\n");
			fclose(bed2bam_fh);
			fclose(Log);
			return 1;
		}
		
// 		printf("Limit buff = %d\n",fin);
		int i=0, j=0;
// 		for (i=0;i<size;i++){
// 			printf("C1 = %d\tC2 = %d\n",Line_coord_1[i],Line_coord_2[i]);
// 		}
		
		
		int limit = size;
		if(!str_eq(Chr_last,Chr_read)){
			limit--;
		}
		Watchdog = count_score(Line_coord_1, Line_coord_2,Count,limit,inicio,fin,cover_score);
		if(Watchdog){
			fprintf(Log,"Failled scoring\n");
			fclose(bed2bam_fh);
			fclose(Log);
			return 1;
		}
// 		printf("adding score = %d\n",cover_score);
// 		for (i=0;i<size;i++){
// 			printf("C1 = %d\tC2 = %d\t Count = %d \n",Line_coord_1[i],Line_coord_2[i],Count[i]);
// 		}
		Watchdog = print_buffer(Line_coord_1, Line_coord_2,Count, limit,&size,Chr_read, fin);
		if(Watchdog){
			fprintf(Log,"Failled printing\n");
			fclose(bed2bam_fh);
			fclose(Log);
			return 1;
		}
		
		
	}
	
// 	printf("size = %d\n",size);
	if(feof(bed2bam_fh)){
		return 0;
	}
	fprintf(Log,"Problem in core function \n");
	return 0;
}

int print_buffer(int *Line_coord_1, int *Line_coord_2,int *Count,int limit,int *size,char *Chr_read,int fin){
	int temp_lim=limit;
	int i=0;
	if(limit < *size){
		for(i=0;i<temp_lim;i++){
			int LC1 = Line_coord_1[i];
			int LC2 = Line_coord_2[i];
			int Count_temp = Count[i];
			printf("%s\t%d\t%d\t%d\n",Chr_read,LC1,LC2,Count_temp);
			recorrer_array(Line_coord_1,*size,i);
			recorrer_array(Line_coord_2,*size,i);
			recorrer_array(Count,*size,i);
			temp_lim--;
			(*size)--;
			i--;
			if(temp_lim<0){
				fprintf(Log,"Error in print_buffer\n");
				return 1;
			}
		}
	}else{
		for(i=0;i<temp_lim;i++){
			int LC1 = Line_coord_1[i];
			int LC2 = Line_coord_2[i];
			int Count_temp = Count[i];
			if(LC2 <= fin){
				printf("%s\t%d\t%d\t%d\n",Chr_read,LC1,LC2,Count_temp);
				recorrer_array(Line_coord_1,*size,i);
				recorrer_array(Line_coord_2,*size,i);
				recorrer_array(Count,*size,i);
				temp_lim--;
				(*size)--;
				i--;
				if(temp_lim<0){
					fprintf(Log,"Error in print_buffer\n");
					return 1;
				}
			}
		}
	}
	return 0;
}

int count_score(int *Line_coord_1, int *Line_coord_2,int *Count,int limit,int inicio,int fin,int cover_score){
	int i = 0;
	int low=0;
	int high=0;
	if(limit >= MAX_ARR){
		fprintf(Log,"Count_limit out of boundaries\n");
		return 1;
	}
	for (i=0; i< limit; i++){
		int LC1 = Line_coord_1[i];
		int LC2 = Line_coord_2[i];
		
		if(LC1 > LC2){
			fprintf(Log,"Coords Error\n");
			return 1;
		}
		
		if(LC1 < fin && LC2 > inicio){
			low = max(LC1,inicio);
			high = min(LC2,fin);
			Count[i] += cover_score*(high-low);
		}
	}
	
	return 0;
}

int min(int a, int b){
	if(a<b){
		return a;
	}
	return b;
}

int max(int a, int b){
	if(a>b){
		return a;
	}
	return b;
}

int fill_buffer(char *Chr_read, char *Chr_last,int fin,int *Line_coord_1, int *Line_coord_2,int *size){
	
	char Chr[Str_len]={0};
	int Coord1=0;
	int Coord2=0;
	if(*size ==0 ){
		if (!feof(bed2bam_fh)){
			if(fscanf(bed2bam_fh,"%s\t%u\t%u\n",Chr,&Coord1,&Coord2)!= 3){
				return 1;
			}
			Line_coord_1[*size]=Coord1;
			Line_coord_2[*size]=Coord2;
			(*size)++;
			str_copy(Chr_read,Chr);
		}
		
	}else{
		Coord1 = Line_coord_1[*size-1];
		Coord2 = Line_coord_2[*size-1];
		str_copy(Chr_last,Chr);
	}
	while(Coord1 < fin && str_eq(Chr,Chr_read) && !feof(bed2bam_fh) ){
		
		
		
		if(fscanf(bed2bam_fh,"%s\t%u\t%u\n",Chr,&Coord1,&Coord2)!= 3){
			return 1;
		}
		Line_coord_1[*size]=Coord1;
		Line_coord_2[*size]=Coord2;
		(*size)++;
		
		if(*size >= MAX_ARR){
			fprintf(Log,"Overload line list\n");
			return 1;
		}
		
	}
	str_copy(Chr,Chr_last);
	
	return 0;
}

int str_copy(char *Str1,char *Str2){
	str_clear(Str2);
	int i = 0;
	while(Str1[i]!=0 && i < Str_len){
		Str2[i]=Str1[i];
		i++;
	}
	if (i >= Str_len){
		fprintf(Log,"String Overload");
		return 1;
	}
	return 0;
}

int str_clear(char *Str){
	int i = 0;
	for(i=0; i < Str_len;i++){
		Str[i]=0;
		i++;
	}
	return 0;
}

int str_eq(char *Str1,char *Str2){
	int i = 0;
	while( i < Str_len){
		if(Str1[i]!=Str2[i]){
			return 0;
		}
		i++;
	}
	return 1;
}

int recorrer_array(int *array, int size, int n){
	
	if(n >= size || n>= MAX_ARR ){
		fprintf(Log,"Error in recorrer_array\n");
		return 1;
	}
	int i=0;
	int lim = size-1;
	for(i=n;i< lim;i++){
		array[i]=array[i+1];
	}
	array[lim]=0;
	return 0;
}
