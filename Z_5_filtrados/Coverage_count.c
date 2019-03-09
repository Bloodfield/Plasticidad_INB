#include <stdio.h>

#define MAX_ARR 500000
#define Log_name "Log.txt"
#define Str_len	40
#define n_chr	460

/*
 * 	Coverage_Count
 * 
 * 	Realiza el conteo de covertura total dentro de un intervalo definido po un bed
 * 
 * 	Input :
 * 		bed2bam 	= Nombre de un archivo BED de coordenadas ordenadas por cromosoma y coordenadas
 * 		stdin		= salida de bedtools genomecov para un genoma en formato bedgraph
 * 	
 * 	Output :
 * 		stdout 		= Salida de las coordenadas de bed2bam con la cuenta de cobertura
 */

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
char Chr_names[n_chr][Str_len]={0};
int Chr_ID[n_chr]={0};
int Chr_len=0;
int visited_chr[n_chr]={0};


//	MAIN
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
	
	//	Algoritmo
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



//	FUNCIONES

//		Sistema

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

// 		Algoritmo

int Coverage_Count(){
	
	//	Variables de estado
	int Watchdog = 0;
	char Chr_last[Str_len]={0};
	char Chr_read[Str_len]={0};
	int inicio = 0;
	int fin = 0;
	int cover_score = 0 ;
	int end_chr = 1;
	
	//	Crea indice de archivo
	Watchdog = create_index();
	if(Watchdog){
		fprintf(Log,"Index failed\n");
		fclose(bed2bam_fh);
		fclose(Log);
		return 1;
	}
	
	//	ALgoritmo
	while( !feof(stdin)){

		//	Get line from stdin
		str_clear(Chr_read);
		if(scanf("%s\t%d\t%d\t%d\n",Chr_read,&inicio,&fin,&cover_score)!= 4){
			fprintf(Log,"Archivo : Not a BED graph format in the stdin file\n");
			return 1;
		}
		
		//	Llena el buffer de los datos de bed2bam
		Watchdog = fill_buffer(Chr_read,Chr_last,fin, &end_chr);
		if(Watchdog){
			fprintf(Log,"Failed to fill buffer\n");
			fclose(bed2bam_fh);
			fclose(Log);
			return 1;
		}
		
		
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
	
	//	Imprime los cromosomas restantes
	int i = 0;
	int counter = 0;
	for (i=0;i < Chr_len; i++){
		if (! visited_chr[i]){
			char Chr[Str_len]={0};
			int Coord1=0;
			int Coord2=0;
			fseek(bed2bam_fh,Chr_ID[i],SEEK_SET);
			if(fscanf(bed2bam_fh,"%s\t%d\t%d\n",Chr,&Coord1,&Coord2)!= 3){
				fprintf(Log,"Err : Coverage_Count : problem reading input\n");
				return 1;
			}
			while( str_eq(Chr,Chr_names[i]) && !feof(bed2bam_fh) ){
				printf("%s\t%d\t%d\t%d\n",Chr,Coord1,Coord2,0);
				if(fscanf(bed2bam_fh,"%s\t%d\t%d\n",Chr,&Coord1,&Coord2)!= 3){
					fprintf(Log,"Err : Coverage_Count : problem reading input\n");
					return 1;
				}
			}
			if( str_eq(Chr,Chr_names[i]) && feof(bed2bam_fh) ){
				printf("%s\t%d\t%d\t%d\n",Chr,Coord1,Coord2,0);
			}
			counter++;
		}
		
	}
	if (counter  > 0){
		fprintf(Log,"Warning : Coverage_Count : 3_0 : %d chromosomes missing in the bed graph input \n",counter);
	}
	return 0;
}

int fill_buffer(char *Chr_read,char *Chr_last, int fin, int *end_chr){
	
	char Chr[Str_len]={0};
	int temp_line = 0;
	int Coord1=0;
	int Coord2=0;
	
	//	En caso de cambiar de cromosoma
	if (! str_eq(Chr_last,Chr_read)){
		
		//	Añade el resto del cromosoma al buffer
		str_copy(Chr_last,Chr);
		while( str_eq(Chr,Chr_last) && !feof(bed2bam_fh) ){
			temp_line=ftell(bed2bam_fh);
			str_clear(Chr);
			if(fscanf(bed2bam_fh,"%s\t%d\t%d\n",Chr,&Coord1,&Coord2)!= 3){
				fprintf(Log,"Err : Fill buffer : problem in line %d\n",temp_line);
				return 1;
			}
			Line_coord_1[size]=Coord1;
			Line_coord_2[size]=Coord2;
			size++;
			
			if(size >= MAX_ARR){
				fprintf(Log,"Overload line list\n");
				return 1;
			}
			
		}
		if (!feof(bed2bam_fh) && size > 0){
			recorrer_tabla(size-1);
		}
		//	Imprime lo que se tiene en el buffer
		
		while(size > 0){
			int LC1 = Line_coord_1[0];
			int LC2 = Line_coord_2[0];
			int Count_temp = Count[0];
			printf("%s\t%d\t%d\t%d\n",Chr_last,LC1,LC2,Count_temp);
			recorrer_tabla(0);
		}
		
		//	Identifica la linea de cromosoma
		
		int idx = search(Chr_read,Chr_names,Chr_len);
		if (idx == -1){
			return 0;
		}
		fseek(bed2bam_fh,Chr_ID[idx],SEEK_SET);
		visited_chr[idx]=1;
		
		//	Resetea Variables
		str_copy(Chr_read,Chr_last);
		*end_chr = 0;
	}
	
	//	Llena los datos faltantes en el buffer dentro del cromosoma
	if ( *end_chr){
		return 0;
	}
	Coord1=0;
	Coord2=0;
	str_clear(Chr);
	str_copy(Chr_read,Chr);
	temp_line = 0;
	
	//		Añade elementos dentro del intervalo
	while(Coord1 < fin && str_eq(Chr,Chr_read) && !feof(bed2bam_fh) ){
		
		temp_line=ftell(bed2bam_fh);
		str_clear(Chr);
		if(fscanf(bed2bam_fh,"%s\t%d\t%d\n",Chr,&Coord1,&Coord2)!= 3){
			fprintf(Log,"Err : Fill buffer : problem in line %d\n",temp_line);
			return 1;
		}
		Line_coord_1[size]=Coord1;
		Line_coord_2[size]=Coord2;
		size++;
		
		if(size >= MAX_ARR){
			fprintf(Log,"Overload line list\n");
			return 1;
		}
		
	}
	
	//		En caso de añadir uno distinto, se elimina
	if (!feof(bed2bam_fh) && size > 0){
		recorrer_tabla(size-1);
		fseek(bed2bam_fh, temp_line,SEEK_SET);
	}
	//		Si se llega al final del cromosoma, levanta la bandera 
	if(feof(bed2bam_fh) || !str_eq(Chr,Chr_read) ){
		*end_chr = 1;
		
	}
	
	return 0;
}

int create_index(){
	
	//	Variables de estado
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
	
	//	suma la cbertura por cada una de las bases entre inicio y fin que intersectan con el buffer
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

//		Basic numerics

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

//		Strings

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

int search(char *Str,char Array[][Str_len],int Array_len){
	int idx=0;
	for(idx=0;idx < Array_len;idx++){
		if(str_eq(Str,Array[idx])){
			return idx;
		}
	}
	return -1;
	
}

//		Arrays

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
