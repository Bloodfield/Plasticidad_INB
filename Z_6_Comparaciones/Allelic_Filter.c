//	RMFilter in C

/*	TODO
 * 	Hacer recvisión de flancos derecho e izquierdo, que en efecto sea izquierdo y derecho
 * 	Hacer guía
 * 	Funcionalizal main
 */
#include <stdio.h>
#include <stdlib.h>
//	0_Core
#define Log_name "Log.txt"
#define Str_len	40
#define n_lines	120000
#define n_chr	460

//	funciones

//	Sistema
int print_help(FILE *Log);

//	Algoritmo
int read_bedcov(FILE *bedcov_fh, char Chr[][Str_len], int *flanco_A, int *flanco_B, int *cuenta, char Chr_names[][Str_len], int *Chr_ID, int *table_len, int *Chr_len, FILE *Log);
int Filter(FILE *Overlap_fh, char Chr[][Str_len], int *flanco_A, int *flanco_B, int *cuenta, char Chr_names[][Str_len], int *Chr_ID, int table_len, int Chr_len,int th, FILE *Log);
int get_cover(char *Chr_search,int inicio,int fin, char Chr[][Str_len], int *flanco_A, int *flanco_B, int *cuenta, char Chr_names[][Str_len], int *Chr_ID, int table_len, int Chr_len);

//	Arrays
int search(char *Str,char Array[][Str_len],int Array_len);
int Array_eq(char *Str1,char *Str2);
int search_from(int query, int *Array,int Array_len, int from);
int copy_Array(char *Str1,char *Str2);

//	Z_Work

int main(int argc, char *argv[]){
	
	//	Head
	
	FILE *Log;
	Log = fopen(Log_name,"a");
	char *bedcov;
	char *Overlap;
	float th=0.0;
		
	if (argc != 4){
		print_help(Log);
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
	
	//	Abrir y probar archivos
	
	FILE *Overlap_fh;
	FILE *bedcov_fh;
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
	
	
	int Whatchdog=0;
	
	//	Leer archivo bedcov
	int flanco_A[n_lines];
	int flanco_B[n_lines];
	int cuenta[n_lines];
	char Chr[n_lines][Str_len];
	
	char Chr_names[n_chr][Str_len];
	int Chr_ID[n_chr];
	
	int table_len=0;
	int Chr_len=0;
	
	Whatchdog = read_bedcov(bedcov_fh, Chr, flanco_A, flanco_B, cuenta, Chr_names, Chr_ID, &table_len, &Chr_len, Log);
	fclose(bedcov_fh);
	
	int i=0;
// 	for (i=0;i<table_len;i++){
// 		printf("Chr = %s\tFlanco A = %d \tFlanco B = %d \tCuenta = %d\n",Chr[i],flanco_A[i],flanco_B[i],cuenta[i]);
// 	}
// 	for (i=0;i<Chr_len;i++){
// 		printf("Chr = %s\tID = %d \n",Chr_names[i],Chr_ID[i]);
// 	}
	
	if (Whatchdog){
		fclose(Overlap_fh);
		fclose(Log);
		return 1;
	}
	
	//	Ejecuta el filtro
	Whatchdog =Filter(Overlap_fh, Chr, flanco_A, flanco_B, cuenta,  Chr_names, Chr_ID, table_len, Chr_len, th, Log);
	
	if (Whatchdog){
		fclose(Overlap_fh);
		fclose(Log);
		return 1;
	}
	
	//	Cierre y salida
	
	fclose(Overlap_fh);
	fclose(Log);
	
	return 0;
}

int print_help(FILE *Log){

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
	fprintf(Log,"$ Allelic_Filter coverarge_del.bed coverarge_left coverarge_right 0.025 \n");
	fprintf(Log,"\n");
	return 0;
}

int read_bedcov(FILE *bedcov_fh, char Chr[][Str_len], int *flanco_A, int *flanco_B, int *cuenta, char Chr_names[][Str_len], int *Chr_ID, int *table_len, int *Chr_len, FILE *Log){
	int i=0,j=0;
	int inicio=0, fin=0, cover_score=0;
	char Chr_previous[Str_len]={0};
	while(!feof(bedcov_fh)){
		
		//	Get line
		char Chr_read[Str_len]={0};
		if(fscanf(bedcov_fh,"%s\t%u\t%u\t%u\n",Chr_read,&inicio,&fin,&cover_score)!= 4){
			fprintf(Log,"Archivo : Not a BED_Coverage format in the bedcov file\n");
			return 1;
		}
		
		copy_Array(Chr_read,Chr[i]);
		flanco_A[i]=inicio;
		flanco_B[i]=fin;
		cuenta[i]=cover_score;
		
		if(! Array_eq(Chr_read,Chr_previous)){
			copy_Array(Chr_read,Chr_names[j]);
			copy_Array(Chr_read,Chr_previous);
			Chr_ID[j]=i;
			j++;
			if (j>= n_chr){
				fprintf(Log,"Archivo : Chr_names_overload\n");
				return 1;
			}
		}
		i++;
		if (j>= n_lines){
			fprintf(Log,"Archivo : n_lines_overload\n");
			return 1;
		}
		
	}
	*table_len=i;
	*Chr_len=j;
	
	return 0;
}

int copy_Array(char *Str1,char *Str2){
	int i=0;
	for(i=0; i< Str_len;i++){
		Str2[i]=Str1[i];
	}
}

int Filter(FILE *Overlap_fh, char Chr[][Str_len], int *flanco_A, int *flanco_B, int *cuenta, char Chr_names[][Str_len], int *Chr_ID, int table_len, int Chr_len,int th, FILE *Log){
	
	float landscape_coverage=0.0, del_coverage=0.0;
	float ratio=0.0;
	int inicio=0, fin=0;
	int overlap_score_del=0;
	int lenght=0;
	int cover_score_del=0,cover_score_left=0, cover_score_right=0;
	
	while(!feof(Overlap_fh)){
		
		//	Get line
		char Chr_search[Str_len]={0};
		if(fscanf(Overlap_fh,"%s\t%u\t%u\t%u\n",Chr_search,&inicio,&fin,&overlap_score_del)!= 4){
			fprintf(Log,"Archivo : Not a BED_Overlap format in the Overlap file\n");
			return 1;
		}
		
		//	search coverage scorings
		
		cover_score_del = get_cover(Chr_search,inicio,fin,Chr,flanco_A,flanco_B,cuenta,Chr_names,Chr_ID,table_len,Chr_len);
		cover_score_left = get_cover(Chr_search,inicio-20,inicio,Chr,flanco_A,flanco_B,cuenta,Chr_names,Chr_ID,table_len,Chr_len);
		cover_score_right = get_cover(Chr_search,fin,fin+20,Chr,flanco_A,flanco_B,cuenta,Chr_names,Chr_ID,table_len,Chr_len);
		
		if (cover_score_del==-1 || cover_score_left==-1 || cover_score_right==-1){
			fprintf(Log,"Cover file not complete\n");
			return 1;
		}
		
		//	compute fraction
		lenght = fin -inicio;
		// 		printf("lenght = %u\n",lenght);
		landscape_coverage=(cover_score_left+cover_score_right)/40.0;
		// 		printf("land cov = %0.6f\n",landscape_coverage);
		del_coverage=cover_score_del/(1.0*lenght);
		// 		printf("del cov = %0.6f\n",del_coverage);
		ratio=del_coverage/landscape_coverage;
		// 		printf("ratio = %0.6f\n",ratio);
		
		//	check and print
		if((ratio < 1-th && ratio >0.5+th) || (ratio <0.5-th && ratio >th)){
			printf("%s\t%u\t%u\t%.5f\t%u\n",Chr_search,inicio,fin,ratio,overlap_score_del);
		}
	}
	
	return 0;
}

int get_cover(char *Chr_search,int inicio,int fin, char Chr[][Str_len], int *flanco_A, int *flanco_B, int *cuenta, char Chr_names[][Str_len], int *Chr_ID, int table_len, int Chr_len){
	int idx=0;
	int Chr_end=0;
	idx = search(Chr_search,Chr_names,Chr_len);
	if(idx <0){return -1;}
	if (idx < Chr_len-1){
		Chr_end=Chr_ID[idx+1];
	}else{
		Chr_end=table_len;
	}
	idx = Chr_ID[idx];
	idx = search_from(inicio,flanco_A,Chr_end,idx);
	if(idx <0){return -1;}
	Chr_end=idx;
	while(flanco_A[Chr_end]==flanco_A[idx] && Chr_end < table_len -1){Chr_end++;}
	if(flanco_A[Chr_end]==flanco_A[idx] && Chr_end==table_len -1){Chr_end++;}
	idx = search_from(fin,flanco_B,Chr_end,idx);
	if(idx <0){return -1;}
	
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
