#include <stdio.h>

#define MAX_ARR 500000
#define Log_name "Log.txt"
#define Str_len	40
#define n_chr	460

//	Sistema
int print_help();

//	Algoritmo
int create_index();
int Coverage_Count();
int fill_buffer(char *Chr_read,char *Chr_last, int inicio, int *end_chr);
int count_score(int inicio,int fin,int cover_score);
int print_buffer(char *Chr_read,int fin);

//	Strings
int str_copy(char *Str1,char *Str2);
int str_clear(char *Str);
int str_eq(char *Str1,char *Str2);
int search(char *Str,char Array[][Str_len],int Array_len);

//	Basic numerics
int min(int a, int b);
int max(int a, int b);

//	Arrays

int recorrer_tabla(int n);

//	Globales

//	*	Archivos
FILE *Log;
FILE *bed2bam_fh;

//	*	Objeto de tabla
int Line_coord_1[MAX_ARR]={0};
int Line_coord_2[MAX_ARR]={0};
int Count[MAX_ARR]={0};
int size = 0;

//	*	Objeto de idexado del archivo bedcov_fh
char Chr_names[n_chr][Str_len];
int Chr_ID[n_chr];
int Chr_len=0;

int main(int argc, char *argv[]){
	Log = fopen(Log_name,"a");
	//	Input 
	
	if (argc != 2){
		fprintf(Log,"argc = %d\n",argc);
		print_help();
		fclose(Log);
		return 1;
	}
	
	fprintf(Log,"Name_bed2bam : %s \n",argv[1]);
	char *bed2bam;
	bed2bam=argv[1];
	
	//	Abrir y probar archivos
	
	
	bed2bam_fh= fopen(bed2bam,"r");
	
	if(!bed2bam_fh ){
		fprintf(Log,"Archivo : %s no existe\n",bed2bam);
		fclose(bed2bam_fh);
		fclose(Log);
		return 1;
	}
	
	//	Programa
	int Watchdog =0;
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
	
	//	Variables de estado
	int Watchdog = 0;
	char Chr_last[Str_len]={0};
	char Chr_read[Str_len]={0};
	int inicio = 0;
	int fin = 0;
	int cover_score = 0 ;
	int end_chr = 1;
	
	Watchdog = create_index();
	if(Watchdog){
		fprintf(Log,"Index failed\n");
		fclose(bed2bam_fh);
		fclose(Log);
		return 1;
	}
	
// 	printf("Echo : Coverage_Count : Chr_ID = [");
// 	int db_i=0;
// 	for(db_i=0;db_i < Chr_len; db_i++){
// 		printf("%d, ",Chr_ID[db_i]);
// 	}
// 	printf("]\n");
	
// 	fprintf(Log,"Coverage_Count : %d chromosomes read\n",Chr_len);
	
	while( !feof(stdin)){
// 		fprintf(Log,"Echo : Coverage_Count : p2_1\n");
		//	Get line from stdin
		if(scanf("%s\t%d\t%d\t%d\n",Chr_read,&inicio,&fin,&cover_score)!= 4){
			fprintf(Log,"Archivo : Not a BED graph format in the stdin file\n");
			return 1;
		}
// 		printf("Echo : Coverage_Count : Line : %s\t%d\t%d\t%d\n",Chr_read,inicio,fin,cover_score);
		
		//	Llena el buffer de los datos de bed2bam
		Watchdog = fill_buffer(Chr_read,Chr_last,fin, &end_chr);
		if(Watchdog){
			fprintf(Log,"Failed to fill buffer\n");
			fclose(bed2bam_fh);
			fclose(Log);
			return 1;
		}
		
// 		printf("Echo : Coverage_Count : p2_1 : Line Coord A = [");
// 		int db_i=0;
// 		for(db_i=0;db_i < size; db_i++){
// 			printf("%d, ",Line_coord_1[db_i]);
// 		}
// 		printf("]\n");
		
		//	Realiza la cuenta de los scores
		Watchdog = count_score(inicio,fin,cover_score);
		if(Watchdog){
			fprintf(Log,"Failled scoring\n");
			fclose(bed2bam_fh);
			fclose(Log);
			return 1;
		}
		
		//	Imprime los valores terminados y purga la tabla
		Watchdog = print_buffer(Chr_read, fin);
		if(Watchdog){
			fprintf(Log,"Failled printing\n");
			fclose(bed2bam_fh);
			fclose(Log);
			return 1;
		}
		
		
	}
	return 0;
}

int fill_buffer(char *Chr_read,char *Chr_last, int fin, int *end_chr){
	
	//	En caso de cambiar de cromosoma
	if (! str_eq(Chr_last,Chr_read)){
		
		//	Print all buffer
		
		while(size > 0){
			int LC1 = Line_coord_1[0];
			int LC2 = Line_coord_2[0];
			int Count_temp = Count[0];
			printf("%s\t%d\t%d\t%d\n",Chr_read,LC1,LC2,Count_temp);
			recorrer_tabla(0);
		}
		
		//	encontrar linea de cromosoma
		
		int idx = search(Chr_read,Chr_names,Chr_len);
		if (idx == -1){
// 			fprintf(Log,"Chr \" %s \" not found in bed2bam\n",Chr_read);
			return 0;
		}
		fseek(bed2bam_fh,Chr_ID[idx],SEEK_SET);
		
		//	Reset Variables
// 		printf("Echo : fill_buffer : p1_1 : %s > %s\n",Chr_read,Chr_last);
// 		printf("Echo : fill_buffer : p1_1 : %s > %s\n",Chr_read,Chr_last);
		str_copy(Chr_read,Chr_last);
		*end_chr = 0;
	}
	
	//	LLena datos faltantes
	char Chr[Str_len]={0};
	int Coord1=0;
	int Coord2=0;
	str_copy(Chr_read,Chr);
	int temp_line = 0;
// 	if(size ==0 ){
// 		if (!feof(bed2bam_fh)){
// 			if(fscanf(bed2bam_fh,"%s\t%d\t%d\n",Chr,&Coord1,&Coord2)!= 3){
// 				return 1;
// 			}
// 			Line_coord_1[size]=Coord1;
// 			Line_coord_2[size]=Coord2;
// 			size++;
// 			// 			str_copy(Chr_read,Chr);
// 		}
// 		
// 	}else{
// 		Coord1 = Line_coord_1[size-1];
// 		Coord2 = Line_coord_2[size-1];
// 		str_copy(Chr_last,Chr);
// 	}
	
	if (! (*end_chr)){
// 		printf("Echo : fill_buffer : p2_1\n");
		while(Coord1 < fin && str_eq(Chr,Chr_read) && !feof(bed2bam_fh) ){
			
			temp_line=ftell(bed2bam_fh);
			if(fscanf(bed2bam_fh,"%s\t%d\t%d\n",Chr,&Coord1,&Coord2)!= 3){
				fprintf(Log,"Err : Fill buffer : problem in line %d\n",temp_line);
				return 1;
			}
			Line_coord_1[size]=Coord1;
			Line_coord_2[size]=Coord2;
			size++;
// 			printf("Echo : fill_buffer : p2_2 : size = %d\n",size);
			
			if(size >= MAX_ARR){
				fprintf(Log,"Overload line list\n");
				return 1;
			}
			
		}
		if (!feof(bed2bam_fh) && size > 0){
			recorrer_tabla(size-1);
			fseek(bed2bam_fh, temp_line,SEEK_SET);
		}
		if(feof(bed2bam_fh) || !str_eq(Chr,Chr_read) ){
			*end_chr = 1;
			
		}
	}
	
	return 0;
}

int create_index(){
	int i=0;
	int inicio=0, fin=0;
	char Chr_previous[Str_len]={0};
	
	while(!feof(bed2bam_fh)){
		
		//	temporal ftell
		int temp=ftell(bed2bam_fh);
		
		//	Get line
		char Chr_read[Str_len]={0};
		if(fscanf(bed2bam_fh,"%s\t%d\t%d\n",Chr_read,&inicio,&fin)!= 3){
			fprintf(Log,"Archivo : Not a BED_Coverage format in the bedcov file\n");
			return 1;
		}
		
		//	si cambió el cromsoma, entonces guardar el elemento de índice
		if(! str_eq(Chr_read,Chr_previous)){
			str_copy(Chr_read,Chr_names[i]);
			str_copy(Chr_read,Chr_previous);
			Chr_ID[i]=temp;
			i++;
			if (i>= n_chr){
				fprintf(Log,"Archivo : Chr_names_overload\n");
				return 1;
			}
		}
	}
	Chr_len=i;
	return 0;
	
}

int print_buffer(char *Chr_read,int fin){
// 	printf("Echo : print_buffer : p1 Size = %d\n",size);
	int i=0;
	for(i=0;i<size;i++){
		int LC1 = Line_coord_1[i];
		int LC2 = Line_coord_2[i];
		int Count_temp = Count[i];
		if(LC2 <= fin){
			printf("%s\t%d\t%d\t%d\n",Chr_read,LC1,LC2,Count_temp);
			recorrer_tabla(i);
			i--;
		}
	}
	return 0;
}

int count_score(int inicio,int fin,int cover_score){
	int i = 0;
	int low=0;
	int high=0;
	if(size >= MAX_ARR){
		fprintf(Log,"Count_limit out of boundaries\n");
		return 1;
	}
	for (i=0; i< size; i++){
		int LC1 = Line_coord_1[i];
		int LC2 = Line_coord_2[i];
		
		if(LC1 >= LC2){
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



int str_copy(char *Str1,char *Str2){
	str_clear(Str2);
	int i = 0;
	for (i=0;i < Str_len;i++){
		Str2[i]=Str1[i];
	}
	return 0;
}

int str_clear(char *Str){
	int i = 0;
	for(i=0; i < Str_len;i++){
		Str[i]=0;
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

int recorrer_tabla( int n){
	
	if(n >= size || n>= MAX_ARR || size <= 0 ){
		fprintf(Log,"Error in recorrer_tabla\n");
		return 1;
	}
	int i=0;
	int lim = size-1;
	for(i=n;i< lim;i++){
		Line_coord_1[i]=Line_coord_1[i+1];
		Line_coord_2[i]=Line_coord_2[i+1];
		Count[i]=Count[i+1];
	}
	Line_coord_1[lim]=0;
	Line_coord_2[lim]=0;
	Count[lim]=0;
	size --;
	return 0;
}

int search(char *Str,char Array[][Str_len],int Array_len){
	int idx=0;
	for(idx=0;idx < Array_len;idx++){
		if(str_eq(Str,Array[idx])){
			return idx;
		}
	}
	return -1;
	
}
