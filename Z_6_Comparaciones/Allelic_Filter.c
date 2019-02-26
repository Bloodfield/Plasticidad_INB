//	RMFilter in C

/*	TODO
 * 	Hacer recvisión de flancos derecho e izquierdo, que en efecto sea izquierdo y derecho
 * 	Hacer guía
 * 	Funcionalizar main
 */

//	0_Core

//	*	Bibliotecas
#include <stdio.h>
#include <stdlib.h>

//	*	Definiciones
#define Log_name "Log.txt"
#define Str_len	40
#define n_lines	100000
#define n_chr	460

//	funciones

//	*	Sistema
int print_help();

//	*	Algoritmo
int index_bedcov( );
int read_bedcov( int *flanco_A, int *flanco_B, int *cuenta, int *table_len, char *Chr_query);
int Filter(  float th);
int get_cover(int inicio,int fin,  int *flanco_A, int *flanco_B, int *cuenta, int table_len);

//	*	Arrays
int search(char *Str,char Array[][Str_len],int Array_len);
int Array_eq(char *Str1,char *Str2);
int search_from(int query, int *Array,int Array_len, int from);
int copy_Array(char *Str1,char *Str2);
int Array_clear(int *Array);

//	global variables

//	*	Objeto de idexado del archivo bedcov_fh
char Chr_names[n_chr][Str_len];
int Chr_ID[n_chr];
int Chr_len=0;

//	*	Archivos
FILE *Log;
FILE *Overlap_fh;
FILE *bedcov_fh;

//	Z_Work

int main(int argc, char *argv[]){
	
	
	Log = fopen(Log_name,"a");
	char *bedcov;
	char *Overlap;
	
	//	Input
	
	float th=0.0;
		
	if (argc != 4){
		print_help();
		fclose(Log);
		return 1;
	}
	
	fprintf(Log,"Name_Overlap : %s \n",argv[1]);
	fprintf(Log,"Name_bedcov : %s \n",argv[2]);
	fprintf(Log,"Threshold : %f \n",atof(argv[3]));
// 	fprintf(Log,"Threshold : %f \n",argv[3]);
	Overlap=argv[1];
	bedcov=argv[2];
// 	sscanf(argv[3],"%f",&th);
	th=atof(argv[3]);
	fprintf(Log,"Threshold : %f \n",th);
	
	//	Abrir y probar archivos
	
	
	Overlap_fh= fopen(Overlap,"r");
	bedcov_fh = fopen(bedcov,"r");
	
	int ret =0;
	if(!Overlap_fh ){
		fprintf(Log,"Archivo : %s no existe\n",Overlap);
		ret= 1;
	}
	
	if(!bedcov_fh ){
		fprintf(Log,"Archivo : %s no existe\n",bedcov);
		ret= 1;
	}
	if (ret){
		fclose(Overlap_fh);
		fclose(bedcov_fh);
		fclose(Log);
		return 1;
	}
	
	//	BEGIN
	
	int Whatchdog=0;
	
	//	*	crea el índice de ftell de bedcov_fh
	Whatchdog = index_bedcov( );
	if (Whatchdog){
		fclose(Overlap_fh);
		fclose(bedcov_fh);
		fclose(Log);
		return 1;
	}
	
	//	*	Ejecuta el filtro
	
	Whatchdog =Filter( th);
	if (Whatchdog){
		fclose(Overlap_fh);
		fclose(bedcov_fh);
		fclose(Log);
		return 1;
	}
	
	//	END
	
	fclose(Overlap_fh);
	fclose(bedcov_fh);
	fclose(Log);
	
	return 0;
}

int print_help(){

	fprintf(Log,"Error. Mostrando ayuda: \n");
	fprintf(Log," \n");
	fprintf(Log,"Allelic_Filter se utiliza para computar la razon de deleción local a un intervalo localizado \n");
	fprintf(Log,"Su uso es: \n");
	fprintf(Log," \n");
	fprintf(Log,"$ Allelic_Filter Archivo_con_deleciones.bed Covertura_Izquierda Covertura_Derecha th \n");
	fprintf(Log," \n");
	fprintf(Log," Los archivos de entrada contienen las coordenadas y la covertura total\n");
	fprintf(Log," Estos se obtienen con la herremienta \"samtools bedcov -j\" de samtools.\n");
	fprintf(Log," \n");
	fprintf(Log," Los flancos tienen que ser de 20 nt de longitud\n");
	fprintf(Log," \n");
	fprintf(Log," th es el umbral para detectar los orígenes germinales o allelicos\n");
	fprintf(Log," Es un número flotante que indica la olgura de las razones para asumir que se encuentra dentro del rango de deleción somática o alélica\n");
	fprintf(Log," \n");
	fprintf(Log,"Ejemplo:\n");
	fprintf(Log,"\n");
	fprintf(Log,"$ Allelic_Filter Overlap.bed bedcov.bed 0.025 \n");
	fprintf(Log,"\n");
	return 0;
}

int index_bedcov(){
	int i=0;
	int inicio=0, fin=0, cover_score=0;
	char Chr_previous[Str_len]={0};
	
	while(!feof(bedcov_fh)){
		
		//	temporal ftell
		int temp=ftell(bedcov_fh);
		
		//	Get line
		char Chr_read[Str_len]={0};
		if(fscanf(bedcov_fh,"%s\t%d\t%d\t%d\n",Chr_read,&inicio,&fin,&cover_score)!= 4){
			fprintf(Log,"Archivo : Not a BED_Coverage format in the bedcov file\n");
			return 1;
		}
		
		//	si cambió el cromsoma, entonces guardar el elemento de índice
		if(! Array_eq(Chr_read,Chr_previous)){
			copy_Array(Chr_read,Chr_names[i]);
			copy_Array(Chr_read,Chr_previous);
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

int read_bedcov(  int *flanco_A, int *flanco_B, int *cuenta, int *table_len, char *Chr_query){
	
	//	asignar lectura desde inicio de cromosoma por medio del índice
	int idx = search(Chr_query,Chr_names,Chr_len);
	if (idx < 0){
		fprintf(Log,"Archivo : Not a BED_Coverage format in the bedcov file\n");
		return 1; 
	}
	fseek(bedcov_fh,Chr_ID[idx],SEEK_SET);
	
	//	variables de estado
	int i=0;
	char Chr_previous[Str_len]={0};
	copy_Array(Chr_query,Chr_previous);
	
	while(!feof(bedcov_fh) && Array_eq(Chr_query,Chr_previous)){
		
		//	Get line
		int inicio=0, fin=0, cover_score=0;
// 		char Chr_read[Str_len]={0};
		if(fscanf(bedcov_fh,"%s\t%d\t%d\t%d\n",Chr_previous,&inicio,&fin,&cover_score)!= 4){
			fprintf(Log,"Archivo : Not a BED_Coverage format in the bedcov file\n");
			return 1;
		}
		
		//	Llenar la tabla
		flanco_A[i]=inicio;
		flanco_B[i]=fin;
		cuenta[i]=cover_score;
		i++;
		if (i>= n_lines){
			fprintf(Log,"Archivo : n_lines_overload\n");
			return 1;
		}
		
	}
	
	if(i==0){
		fprintf(Log,"Archivo : Chromosomes not complete\n");
		return 1;
	}
	
	if(!feof(bedcov_fh)){
		flanco_A[i]=0;
		flanco_B[i]=0;
		cuenta[i]=0;
		*table_len=i-1;
		
		return 0;
	}
	
	*table_len=i;
	return 0;
}

int copy_Array(char *Str1,char *Str2){
	int i=0;
	for(i=0; i< Str_len;i++){
		Str2[i]=Str1[i];
	}
}

int Filter(float th){
	
	//	Objeto de tabla de scoring de bedvoc
	int flanco_A[n_lines];
	int flanco_B[n_lines];
	int cuenta[n_lines];
	
	int table_len=0;
	
	//	variables de estado
	char chr_act[Str_len]={0};
	
	while(!feof(Overlap_fh)){
		
		//	Get line
		int inicio=0, fin=0;
		int overlap_score_del=0;
		char chr_new[Str_len]={0};
		if(fscanf(Overlap_fh,"%s\t%d\t%d\t%d\n",chr_new,&inicio,&fin,&overlap_score_del)!= 4){
			fprintf(Log,"Archivo : Not a BED_Overlap format in the Overlap file\n");
			return 1;
		}
		
		//	llena la tabla para el chromosoma actual si es necesario
		if(! Array_eq(chr_new, chr_act)){
			Array_clear(flanco_A);
			Array_clear(flanco_B);
			Array_clear(cuenta);
			table_len=0;
			read_bedcov(flanco_A, flanco_B, cuenta, &table_len,chr_new);
			copy_Array(chr_new, chr_act);
		}
		
		//	search coverage scorings, se tiene garantizado que solo hay el cromosoma actual está en el buffer
		int cover_score_del = get_cover(inicio,fin,flanco_A,flanco_B,cuenta,table_len);
		int cover_score_left = get_cover(inicio-20,inicio,flanco_A,flanco_B,cuenta,table_len);
		int cover_score_right = get_cover(fin,fin+20,flanco_A,flanco_B,cuenta,table_len);
		
		int Whatchdog=(cover_score_del==-1 || cover_score_left==-1 || cover_score_right==-1);
		if (Whatchdog){
			fprintf(Log,"Cover file not complete\n");
			return 1;
		}
		
		//	compute fraction
		int lenght = fin -inicio;
		// 		printf("lenght = %u\n",lenght);
		float landscape_coverage=(cover_score_left+cover_score_right)/40.0;
		// 		printf("land cov = %0.6f\n",landscape_coverage);
		float del_coverage=cover_score_del/(1.0*lenght);
		// 		printf("del cov = %0.6f\n",del_coverage);
		float ratio=del_coverage/landscape_coverage;
		// 		printf("ratio = %0.6f\n",ratio);
		
		//	check and print
		Whatchdog = (ratio < 1.0-th && ratio >0.5+th) || (ratio <0.5-th && ratio > th);
		if(Whatchdog){
			printf("%s\t%d\t%d\t%.5f\t%d\n",chr_new,inicio,fin,ratio,overlap_score_del);
		}
	}
	
	return 0;
}

int get_cover(int inicio,int fin, int *flanco_A, int *flanco_B, int *cuenta, int table_len){
	int idx=0;
	
	//	busca la primera coordenada
	idx = search_from(inicio,flanco_A,table_len,0);
	if(idx <0){
		fprintf(Log,"get_cover (1) : bedcov not completed, %d not found \n",inicio);
		return -1;
	}
	
	//	busca límite de busqueda para 2da coordenada
	int coord_end=idx;
	while(flanco_A[coord_end]==flanco_A[idx] && coord_end < table_len -1){coord_end++;}
	if(flanco_A[coord_end]==flanco_A[idx] && coord_end==table_len -1){coord_end++;}
	
	//	busca la segunda coordenada
	idx = search_from(fin,flanco_B,coord_end,idx);
	if(idx <0){
		fprintf(Log,"get_cover (2) : bedcov not completed, %d not found \n",fin);
		return -1;
	}
	
	//	regresa el valor de la cuenta
	return cuenta[idx];
}

int search(char *Str,char Array[][Str_len],int Array_len){
	int idx=0;
	for(idx=0;idx < Array_len;idx++){
		if(Array_eq(Str,Array[idx])){
			return idx;
		}
	}
	return -1;
	
}

int Array_eq(char *Str1,char *Str2){
	int i=0;
	for (i=0;i<Str_len;i++){
		if(Str1[i]!=Str2[i]){
			return 0;
		}
	}
	return 1;
}

int search_from(int query, int *Array,int Array_len, int from){
	int idx=0;
	for(idx=from;idx < Array_len;idx++){
		if(query==Array[idx]){
			return idx;
		}
	}
	return -1;
}

int Array_clear(int *Array){
	int i=0;
	for(i=0; i<n_lines;i++){
		Array[i]=0;
	}
}
