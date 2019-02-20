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
#include <stdlib.h>

#define Log_name "Log.txt"
#define Array_Size 1200
#define Classes_max_size 900
#define Clases_adj_max 900
#define Adj_list_max_size 850
#define Str_len 300
#define n_chr 460

//	Arrays
int queue(int *array,int element, int end);
int copy_array(int *array_a,int *array_b,int size);
int add_to_order_array(int *array, int *len , int elem);
int array_intersect_self(int *self_array,int *array_len, int *second,int second_len);
int clear_array(int *array, int size);
int search_int(int query, int *Array,int Array_len);
int remove_in_order_array(int *cluster,int *cluster_length,int node);

//	Strings
int Array_eq(char *Str1,char *Str2);
int copy_Str(char *Str1,char *Str2);
int search(char *Str,char Array[][Str_len],int Array_len);
int clear_Chr_print();

//	Calculos de programa
int percentage_overlap(int FAR,int FBR , int FAQ,int FBQ );
int test_cc(int *cluster,int cluster_length);

//	Sistema
int err_message();

//	Funciones de proceso
int Del_Overlap(int overlap_th, int score_th);
int index_file();
int read_line(  int overlap_th,int *Next_X_Index, int *next_line,int FA0,int FB0);
int read_line_2(  int overlap_th,int *Next_X_Index, int *next_line, int FA0,int FB0);
int conected_component(int *cluster,int *cluster_length, int node ,int score_th);
int print_cluster(int *cluster, int cluster_length, int score_th);
int add_element( int line_number, int overlap_th,int FA, int FB);
int complete_data( int overlap_th, int *next_line);
int purge_classes( int ID);

//	Set operations
int in_class(int *cluster,int cluster_length);
int is_subeq(int *Arr1,int size1, int *Arr2, int size2);

//	Funciones basicas
int max(int a , int b);
int min(int a , int b);

//	Listas de adyacencia
int add_edge_adj_list( int nodo_a, int nodo_b);

//	Variables globales
	//	Archivos
FILE *Log;
FILE *In_file;

//	*	Objeto grafo
int flanco_A[Array_Size]={0};
int flanco_B[Array_Size]={0};
int adj_list[Array_Size][Adj_list_max_size]={0};
int adj_list_sizes[Array_Size]={0};
int factor[Array_Size]={0};
int fin=0;
int line_ID[Array_Size]={0}; //	Indice de lectura

//	*	Objeto de clases
int Classes[Classes_max_size][Clases_adj_max]={0};
int Classes_sizes[Clases_adj_max];
int Classes_size=0;

//	*	Objeto de idexado del archivo bedcov_fh
char Chr_names[n_chr][Str_len];
int Chr_ID[n_chr];
int Chr_len=0;

//	*	Nombre del cromosoma leido

char Chr_print[Str_len]={0};

int main(int argc, char *argv[]){
	
	
	Log = fopen(Log_name,"a");
	
	char *In_file_name;
	int score_th   = 1;
	int overlap_th = 1;
	
	//	Entrada de datos
	
	if (argc != 4){
		err_message();
		fclose(Log);
		return 1;
	}
	
	fprintf(Log,"Name : %s \n",argv[1]);
	In_file_name=argv[1];
	sscanf (argv[2],"%u",&score_th);
	sscanf (argv[3],"%u",&overlap_th);
	if(overlap_th>100){
		fprintf(Log,"Porcentaje de overlap no válido \n");
		err_message();
		fclose(Log);
		return 1;
	}
	
	// Ejecución del programa
	In_file= fopen(In_file_name,"r");
	if(!In_file){
		fprintf(Log,"Archivo : %s no existe\n",In_file_name);
		return 1;
	}
	
	//	Indice de cromosomas en la entrada
	int Overflow = index_file();
	if(Overflow ==1){
		fprintf(Log,"Failed to make index\n");
		return 1;
	}
	int i =0;
// 	fprintf(Log,"Chr\tID\n");
// 	for(i=0;i< Chr_len; i++){
// 		fprintf(Log,"%s\t%d\n",Chr_names[i],Chr_ID[i]);
// 	}
	
	Overflow = Del_Overlap(overlap_th,score_th);
	
	// Cierre de archivos
	fclose(In_file);
	fclose(Log);
	
	if (Overflow==1){return 1;}
	
	
	return 0;
}

/*
 * 	PROGRAMAS DE SISTEMA
 * 	
 * 	
 */

//	Arrays

int queue(int *array,int element, int end){
	if(end + 1 >= Array_Size){
		fprintf(Log,"Err queue = Stack Overflow\n");
		return 1;
	}
	
	array[end] = element;
	
	return 0;
}

int copy_array(int *array_a,int *array_b,int size){
	int i=0;
	if(size >= Array_Size){
		fprintf(Log,"Err copy_array = Stack Overflow\n");
		return 1;
	}
	for(i=0;i<size; i++){
		array_b[i]=array_a[i];
	}
	// 	printf("Echo 6 \n");
	return 0;
}

int add_to_order_array(int *array, int *len , int elem){
	int i =0;
	if(*len >= Array_Size){
		fprintf(Log,"Err add to order array = Stack Overflow\n");
		return 1;
	}
	for (i=0; i < *len && array[i]<elem;i++){}
	int temp = elem;
	if(temp!=array[i]){
		for (i; i <= *len ;i++){
			temp ^= array[i];
			array[i] ^= temp;
			temp ^= array[i];
		}
		(*len)++;
	}else if(i >= *len){
		(*len)++;
	}
	// 	printf("Echo 7 \n");
	return 0;
}

int array_intersect_self(int *self_array,int *array_len, int *second,int second_len){
	
	if(*array_len>Adj_list_max_size || second_len >Adj_list_max_size){return 1;}
	int temp[Adj_list_max_size]={0};
	int temp_len=0;
	int i=0,j=0;
	while(i < *array_len && j<second_len){
		if(self_array[i]==second[j]){
			add_to_order_array(temp,&temp_len,self_array[i]);
			i++;
			j++;
		}else if(self_array[i]<second[j]){
			i++;
		}else if(self_array[i]>second[j]){
			j++;
		}
	}
	
	for (i=0; i< *array_len;i++){
		self_array[i]=0;
	}
	(*array_len)=temp_len;
	copy_array(temp,self_array,*array_len);
	// 	printf("Echo 8 \n");
	return 0;
}

int clear_array(int *array, int size){
	int i=0;
	if(size > Array_Size){
		fprintf(Log,"Err clear array = Stack Overflow\n");
		return 1;
	}
	for (i=0;i<size;i++){
		array[i]=0;
	}
	return 1;
}

int search_int(int query, int *Array,int Array_len){
	int idx=0;
	for(idx=0;idx < Array_len;idx++){
		if(query==Array[idx]){
			return idx;
		}
	}
	return -1;
	
}

int remove_in_order_array(int *array, int *len , int elem){
	int i =0;
	int len_1 = (*len) -1;
	if(*len >= Array_Size){
		fprintf(Log,"Err remove order array= Stack Overflow\n");
		return 1;
	}
	for (i=0; i < len_1 && array[i]<elem;i++){}
	if(elem == array[i] && i < *len){
		for (i; i < len_1 ;i++){
			array[i] = array[i+1];
		}
		array[i] = 0;
		(*len)--;
		return 0;
	}else{
		return -1;
	}
}

//	Strings

int Array_eq(char *Str1,char *Str2){
	int i=0;
	for (i=0;i<Str_len;i++){
		if(Str1[i]!=Str2[i]){
			return 0;
		}
	}
	return 1;
}

int search(char *Str,char Array[][Str_len],int Array_len){
	int idx=0;
	if(Array_len >= Array_Size){
		fprintf(Log,"Err search = Stack Overflow\n");
		return 1;
	}
	for(idx=0;idx < Array_len;idx++){
		if(Array_eq(Str,Array[idx])){
			return idx;
		}
	}
	return -1;
	
}

int copy_Str(char *Str1,char *Str2){
	int i=0;
	for(i=0; i< Str_len;i++){
		Str2[i]=Str1[i];
	}
}

int clear_Chr_print(){
	int i =0;
	for (i=0; i < Str_len;i++){
		Chr_print[i]=0;
	}
}

//	Calculos de programa

int test_cc(int *cluster,int cluster_length){
	
	int i=0;
	if(cluster_length >= Array_Size){
		fprintf(Log,"Err test cc = Stack Overflow\n");
		return 1;
	}
	for (i=0;i<cluster_length;i++){
		int node1 = cluster[i];
		int temp[Array_Size]={0};
		int temp_len = cluster_length;
		copy_array(cluster,temp,temp_len);
		remove_in_order_array(temp,&temp_len,node1);
		if(! is_subeq(temp,temp_len,adj_list[i],adj_list_sizes[i])){
			return 0;
		}
	}
	
	return 1;
}

int percentage_overlap(int FAR,int FBR , int FAQ,int FBQ ){
// int percentage_overlap(int reference , int query ){
// 	if(reference >= Array_Size || query >= Array_Size){
// 		fprintf(Log,"Err percentaje calculation = Stack Overflow\n");
// 		return 1;
// 	}
// 	int FAR = flanco_A[reference];
// 	int FBR = flanco_B[reference];
// 	int FAQ = flanco_A[query];
// 	int FBQ = flanco_B[query];
	// 	printf("ref \t %u\t%u\nque\t%u\t%u\n",FAR,FBR,FAQ,FBQ);
	if (FBQ < FAR || FAQ > FBR){return 0;}
	int divisor = FBR - FAR;
	int temp= 0;
	if (FAQ < FAR){
		if(FBQ < FBR){
			temp = (FBQ-FAR );
			// 			printf("\ttemp1=%u\n",temp);
		}else{
			return 100;
		}
	}else{
		if(FBQ < FBR){
			temp = (FBQ-FAQ );
			// 			printf("\ttemp3=%u\n",temp);
		}else{
			temp = (FBR-FAQ );
			// 			printf("\ttemp4=%u\n",temp);
		}
	}
	// 	printf("temp=%u \tdivisor = %u\n",temp,divisor);
	temp *= (100.0/divisor);
	return temp;
}

//	Sistema

int err_message(){
	fprintf(Log,"Se tiene que escribir un archivo bed para analizar \n");
	fprintf(Log,"El bed tiene que ser un bed separado por cromosomas y estar ordenado de menor a mayor en base a la primera coordenada \n");
	fprintf(Log,"$ Del_Overlap [Name.bed] overlap_th score_th\n");
	fprintf(Log,"Donde overlap_th es el porcentaje entero mínimo para decir que dos lecturas se sobrelapan y\n");
	fprintf(Log,"\t score_th es la mínima cantidad de sobrelapes para reportar una deleción\n");
	fprintf(Log,"Ejemplo:\n");
	fprintf(Log,"\n");
	fprintf(Log,"$ Del_Overlap $HOME/Deletions_Chr1.bed 50 5\n");
	fprintf(Log,"\n");
	fprintf(Log,"El output se reporta en standard output \n");
	fprintf(Log,"\n");
	return 0;
}

//	Funciones de proceso

int index_file(){
	int i=0;
	int inicio=0, fin=0;
	char Chr_previous[Str_len]={0};
	
	while(!feof(In_file)){
		
		//	temporal ftell
		int temp=ftell(In_file);
		
		//	Get line
		char Chr_read[Str_len]={0};
		char Dummy2[Str_len]={0};
		char Dummy3[Str_len]={0};
		char Dummy4[Str_len]={0};
		if(fscanf(In_file,"%s\t%d\t%d\t%s %s %s\n",Chr_read,&inicio,&fin,Dummy2,Dummy3,Dummy4)!= 6){
			fprintf(Log,"%s\t%d\t%d\t%s\n",Chr_read,inicio,fin,Dummy2);
			fprintf(Log,"Archivo : Not a BED file as input\n");
			return 1;
		}
		
		
		//	si cambió el cromsoma, entonces guardar el elemento de índice
		if(! Array_eq(Chr_read,Chr_previous)){
			copy_Str(Chr_read,Chr_names[i]);
			copy_Str(Chr_read,Chr_previous);
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

int Del_Overlap(int overlap_th, int score_th){
	
	//	Variables de estado
	int current_line = -1, next_line=0;
	int Overflow = 0;
	
	while( current_line != next_line && Overflow < 2 ){
		
		//	Clear Variables
		clear_array(flanco_A,Array_Size);
		clear_array(flanco_B,Array_Size);
		int i=0;
		for (i=0;i<Array_Size;i++){
			clear_array(adj_list[i],Adj_list_max_size);
		}
		clear_array(adj_list_sizes,Array_Size);
		clear_array(factor,Array_Size);
		clear_array(line_ID,Array_Size);
		fin=0;
		
		//	Coloca la linea de lectura correspondiente
		fseek(In_file,next_line,SEEK_SET);
		current_line= next_line;
		clear_Chr_print();
		//	Completa las lecturas con los datos del archivo
		if(Overflow < 2){
			Overflow=complete_data( overlap_th,&next_line);
// 			printf("Echo 2.1 \n");
		}
		//	Si algo falló, regresa el fallo
		if(Overflow==2){
			fprintf(Log,"Complete data return non zero sataus\n");
			return 1;
		}
		
// 		fprintf(Log,"ID list \tF A \t F B\n");
// 		for(i=0;i< fin;i++){
// 			fprintf(Log,"%d\t%d\t%d\n",line_ID[i],flanco_A[i],flanco_B[i]);
// 		}
		
		//	Determinar el nodo de referencia
		int current_node = search_int(current_line,line_ID,fin) ;
		if(current_node == -1){
			fprintf(Log,"line ID %d not found\n",current_line);
			fprintf(Log,"number of nodes = %d\n",fin);
			return 1;
		}
		
		
// 		printf("Echo 3.0 \n");
		//	Obten los maximales del grafo
		int finished = 0;
		int cluster[Array_Size]={0};
		int cluster_length=0;
		
// 		printf("call connected conected_component of %d nodes \n",fin);
		finished = conected_component(cluster,&cluster_length,current_node, score_th);
		if(finished ==2){
			
			return 1;
		}
		
		//	Purgar Clases
		int test = purge_classes(current_line);
		if(test!=0){
			fprintf(Log,"Purge class failed\n");
			return 1;
			
		}
	}
	
	return 0;
	
}

int purge_classes( int ID){
	
	if(Classes_size >= Array_Size){
		fprintf(Log,"Err purge classes = Stack Overflow\n");
		return 1;
	}
	
	int size = Classes_size;
	int i=0, j=0;
	for(i=0;i<size;i++){
		for (j=0;j<Classes_sizes[i] && Classes[i][j] < ID; j++){}
		// 			printf("Purge j1 = %u \n",j);
		if(Classes[i][j] == ID && Classes_sizes[i]>0 ){
			int lim=Classes_sizes[i]-1;
			for (j;j<lim;j++){
				Classes[i][j] = Classes[i][j+1];
			}
			Classes[i][j] = 0;
			Classes_sizes[i]--;
		}
	}
	//		Eliminar clases vacías
	for(i=0;i<size && Classes_sizes[i]!=0;i++){}	//	Busca una linea vacía
	if(Classes_sizes[i]==0 && i<size){
		int offset=1;
		//	Realiza corrimiento de datos
		while(Classes_sizes[i+offset]==0 && i+offset<size){
			offset++;
		}
		for(i;i<size-offset;i++){
			Classes_sizes[i]=Classes_sizes[i+offset];
			for(j=0;j<Classes_sizes[i];j++){
				Classes[i][j]=Classes[i+offset][j];
			}
			while(Classes_sizes[i+offset]==0 && i+offset< size-1){
				offset++;
			}
			if(Classes_sizes[i+offset]==0 && i+offset==size-1){
				offset++;
			}
		}
		
		//	Vuelve 0 las últimas
		for(i;i<size;i++){
			for (j;j<Classes_sizes[i];j++){
				Classes[i][j]=0;
			}
			Classes_sizes[i]=0;
		}
		
		//	Actualiza Classes_size
		Classes_size -= offset;
	}
	return 0;
}

int complete_data(int overlap_th, int *next_line){
	
	//	Find file line for chromosome
	char Chr_current[Str_len]={0};
	char Dummy2[Str_len]={0};
	char Dummy3[Str_len]={0};
	char Dummy4[Str_len]={0};
	int begin=0, end=0;
	if(fscanf(In_file,"%s\t%d\t%d\t%s %s %s\n",Chr_current,&begin,&end,Dummy2,Dummy3,Dummy4)!= 6){
		fprintf(Log,"Error in next line = %d\n ",*next_line);
		return 1;
	}
	fprintf(Log,"complete data : Line : %s\t%d\t%d\t%s %s %s\n",Chr_current,begin,end,Dummy2,Dummy3,Dummy4);
	int FA0= begin;
	int FB0= end;
	int idx = search(Chr_current,Chr_names,Chr_len);
	if (idx < 0){
		fprintf(Log,"Archivo : chromosome %s not found\n",Chr_current);
		return 1; 
	}
	fseek(In_file,Chr_ID[idx],SEEK_SET);
	copy_Str(Chr_names[idx],Chr_print);
	fprintf(Log,"complete data : Chr_print = %s\n",Chr_print);
	fprintf(Log,"Complete Data : chr line = %ld\n",ftell(In_file));
	
	//	Skip FAN < FAN_min
	int temp_line=0;
	int lim_buffer = end - ((end - begin)*(100.0/overlap_th))-1;
	fprintf(Log,"buffer limit = %d\n",lim_buffer);
	begin=0;
	while(begin < lim_buffer ){
		temp_line = ftell(In_file);
		fprintf(Log,"complete data : cicle in skip\n");
		if(fscanf(In_file,"%s\t%d\t%d\t%s %s %s\n",Chr_current,&begin,&end,Dummy2,Dummy3,Dummy4)!= 6){
			fprintf(Log,"%s\t%d\t%d Did not found it minimum\n ",Chr_current,begin,end);
			return 1;
		}
	}
	
	fprintf(Log,"complete data : last line after FAN_min:\n");
	fprintf(Log,"complete data : %s\t%d\t%d\t%s %s %s\n",Chr_current,begin,end,Dummy2,Dummy3,Dummy4);
	fseek(In_file,temp_line,SEEK_SET);
	fprintf(Log,"Complete Data : first line = %ld\n",ftell(In_file));
	
	//	Add FAN < FA_read
	int Overflow = 0;
	temp_line=-1;
	while(temp_line<(*next_line)  && Overflow == 0){
		Overflow=read_line(overlap_th,&begin,&temp_line,FA0,FB0);
		if (Overflow >1){
			fprintf(Log,"Read : Problem getting line\n");
			return 2;
		}
		//	Revisa el orden
		if((fin)>1){
			if(flanco_A[(fin)-1] < flanco_A[(fin)-2]){
				fprintf(Log,"File : Bed file not sorted\n");
				// 				printf("Bedfile not sorted\n");
				fclose(Log);
				return 2;
			}
			if(flanco_A[(fin)-1] == flanco_A[0] && flanco_B[(fin)-1] == flanco_B[0]){
				fprintf(Log,"Soft: Error in repetaed coords previous (%d,%d)\nSize=%d", flanco_A[0],flanco_B[0],fin);
				fclose(Log);
				return 2;
			}
		}
	}
	
// 	fprintf(Log," %d temp == %d next\n",temp_line,*next_line);
	//	Add FAN < FAN_max	(get next_line)
	
	lim_buffer = FA0 +((FB0-FA0)*(overlap_th/100.0))+1;	// +1 es por el redondeo
	int flag = 1;
	fprintf(Log,"buffer limit = %d\n",lim_buffer);
	
	while(begin < lim_buffer && Overflow == 0){
// 		fprintf(Log," %d begin < lim_buffer %d \n",begin,lim_buffer);
		
		//	Obten la información de la linea
		Overflow=read_line_2(overlap_th,&begin,&temp_line,FA0,FB0);
		if(flag && temp_line!=(*next_line)){	// cambia el valor de next line al primer cambio
			flag=0;
			(*next_line)=temp_line;
// 			fprintf(Log," Echo complete data : change next line\n");
		}
		if (Overflow >1){
			fprintf(Log,"Overlflow in read_lin\n");
			return 2;
		}
		//	Revisa el orden
		if((fin)>1){
			if(flanco_A[(fin)-1] < flanco_A[(fin)-2]){
				fprintf(Log,"Bed file not sorted\n");
				err_message();
				fclose(Log);
				return 2;
			}
		}
	}
	
	return Overflow;
}

int read_line(int overlap_th,int *X_Index, int *next_line, int FA0, int FB0){
	
	char Dummy1[Str_len]={0};
	char Dummy2[Str_len]={0};
	char Dummy3[Str_len]={0};
	char Dummy4[Str_len]={0};
	int FAN=0,FBN=0;
	int Overflow= 0;
// 	int temp_line=ftell(In_file);
	//	Leer linea
	int line_number = ftell(In_file);
	if(fscanf(In_file,"%s\t%d\t%d\t%s %s %s\n",Dummy1,&FAN,&FBN,Dummy2,Dummy3,Dummy4)!= 6){
		
		return 1;
	}
	// 	printf("%s\t%u\t%u\t%s\n",Dummy1,FAN,FBN,Dummy2);
	(*X_Index)=FAN;
	
	if(fin == 0){
		//	Lee primera linea
		Overflow = add_element( line_number, overlap_th,FAN, FBN);
		factor[0]++;
		(*next_line)=line_number;
// 		fprintf(Log,"Add first element : \t%d\t%d\t%d\n",line_number,FAN,FBN);
		// 		printf("Echo 4.1\tFAN=%u\tFBN=%u\tF_AN=%u\tF_BN=%u\n",FAN,FBN,flanco_A[0],flanco_B[0]);
	}else{
		//	Condiciones de addición:
		int fin_1=(fin)-1;
		int FAP=flanco_A[fin_1];
		int FBP=flanco_B[fin_1];
// 		int FA0=flanco_A[0];
// 		int FB0=flanco_B[0];
		
		//	la próxima linea diferente avalor anterior:
		if(FAP !=FAN || FBP !=FBN){
			fprintf(Log," Read line : change next line\n");
			fprintf(Log,"%d\t%s\t%d\t%d\t%s %s %s\n",line_number,Dummy1,FAN,FBN,Dummy2,Dummy3,Dummy4);
			(*next_line)=line_number;
		}
		
		int condition1=0,condition2=0;
		condition1 = (percentage_overlap(FA0,FB0,FAN,FBN) > overlap_th );
		condition2 = (percentage_overlap(FAN,FBN,FA0,FB0) > overlap_th );
// 		condition1 =(((FBN-FA0)*(100.0/overlap_th)) > (FB0-FA0) ); //	restringe el minimo del flanco derecho.
// 		condition2 =(((100.0/overlap_th)*(FB0-FA0))>(FBN-FAN));	//	restringe el máximo del flanco derecho.
		fprintf(Log,"Add element : cond 1\t%d\tcond2\t%d\n",condition1,condition2);
		if(condition1 && condition2){
			
			if(FAP !=FAN || FBP !=FBN){
				Overflow = add_element(line_number, overlap_th,FAN, FBN);
// 				fprintf(Log,"Add element : \t%d\t%d\t%d\n",line_number,FAN,FBN);
				fin_1=(fin)-1;	// porque add_element cambia el valor de fin
				// 			printf("Echo 4.3\n");
			}
			factor[fin_1]++;
		}
		// 		printf("Echo 4.2\n");
	}
	
	if (Overflow > 0){return  2;}
	
	return 0;
}

int read_line_2(int overlap_th,int *X_Index, int *next_line, int FA0,int FB0){
	
	char Chr_line[Str_len]={0};
	char Dummy2[Str_len]={0};
	char Dummy3[Str_len]={0};
	char Dummy4[Str_len]={0};
	int FAN=0,FBN=0;
	int Overflow= 0;
// 	int temp_line=ftell(In_file);
	//	Leer linea
	int line_number = ftell(In_file);
	if(fscanf(In_file,"%s\t%d\t%d\t%s %s %s\n",Chr_line,&FAN,&FBN,Dummy2,Dummy3,Dummy4)!= 6){
		return 1;
	}
	if(! Array_eq(Chr_line,Chr_print)){
		fprintf(Log,"Chr end \n");
		(*next_line)=line_number;
		return 1;
	}
	// 	printf("%s\t%u\t%u\t%s\n",Chr_line,FAN,FBN,Dummy2);
	(*X_Index)=FAN;
	
	if(fin == 0){
		//	Lee primera linea
		Overflow = add_element( line_number, overlap_th,FAN, FBN);
		factor[(fin)-1]++;
		(*next_line)=line_number;
		// 		printf("Echo 4.1\tFAN=%u\tFBN=%u\tF_AN=%u\tF_BN=%u\n",FAN,FBN,flanco_A[0],flanco_B[0]);
	}else{
		//	Condiciones de addición:
		int fin_1=(fin)-1;
		int FAP=flanco_A[fin_1];
		int FBP=flanco_B[fin_1];
// 		int FA0=flanco_A[0];
// 		int FB0=flanco_B[0];
		int FANI=FAN, FBNI=FBN;
		
		//	la próxima linea diferente la enterior:
		if(FAP !=FAN || FBP !=FBN){
			(*next_line)=line_number;
		}
		
		int condition1=0,condition2=0;
		condition1 = (percentage_overlap(FA0,FB0,FAN,FBN) > overlap_th );
		condition2 = (percentage_overlap(FAN,FBN,FA0,FB0) > overlap_th );
// 		condition1 =(((FBNI-FANI)*100) > ((FB0-FA0)*overlap_th) ); //	Evita añadir cosas muy pequeñas
// 		condition2 =((100*(FB0-FANI))>(overlap_th*(FBNI-FANI)));	// evita añadir cosas muy grandes
		if(condition1 && condition2){
			
			if(FAP !=FAN || FBP !=FBN){
				Overflow = add_element(line_number, overlap_th,FAN, FBN);
				fin_1=(fin)-1;
				// 			printf("Echo 4.3\n");
			}
			factor[fin_1]++;
		}
		// 		printf("Echo 4.2\n");
	}
	
	if (Overflow > 0){return  2;}
	
	return 0;
}

int add_element(int line_number ,int overlap_th,int FA, int FB){
	if(queue(flanco_A,FA,fin)){
		fprintf(Log,"Queue error in add element \n");
		return 1;
	}
	if(queue(flanco_B,FB,fin)){
		fprintf(Log,"Queue error in add element \n");
		return 1;
	}
	if(queue(line_ID,line_number,fin)){
		fprintf(Log,"Queue error in add element \n");
		return 1;
	}
	int i = 0;
	for (i=0;i< fin ; i++){
		int per_ida = percentage_overlap(flanco_A[i],flanco_B[i],flanco_A[fin],flanco_B[fin]);
		int per_regreso = percentage_overlap(flanco_A[fin],flanco_B[fin],flanco_A[i],flanco_B[i]);
		
		if(per_ida > overlap_th && per_regreso > overlap_th){
			if(adj_list_sizes[i] >= Adj_list_max_size||adj_list_sizes[fin] >= Adj_list_max_size){
				fprintf(Log,"AjList overflow\n");
				// 				printf("AjList overflow\n");
				return 1;
				
			}
			
			add_edge_adj_list( i,fin);
			// 			matrix[i][fin]=matrix[fin][i]=1;
		}
	}
	fin++;
	return 0;
}

int print_cluster(int *cluster, int cluster_length, int score_th){
	if(cluster_length >= Array_Size){
		fprintf(Log,"Queue error in add element \n");
		return 1;
	}
	int coord_a=flanco_A[0];
	int coord_b=flanco_B[0];
	int i = 0, score = 0;
	// 	printf("cluster _len= %u\n",cluster_length);
	for (i=0;i<cluster_length;i++){
		int index = cluster[i];
		if (index >= fin){return -1;}
		// 		printf("score= %u\n",factor[index]);
		score += factor[index];
		coord_a = max(coord_a,flanco_A[index]);
		coord_b = min(coord_b,flanco_B[index]);
	}
	if(score >= score_th){
		printf("%s\t%d\t%d\t%d\n",Chr_print,coord_a,coord_b,score);
	}
	
	return 0;
}

int conected_component(int *cluster,int *cluster_length, int node ,int score_th){ // visited not needed
	
	//	añadir nodo
	if(add_to_order_array(cluster,cluster_length,node)){
		fprintf(Log,"Err connected conected_component init \n");
		return 1;
	}
	
// 	printf("cluster = [");
// 	int i_db= 0;
// 	for (i_db = 0; i_db < *cluster_length;i_db++){
// 		printf("%d,",cluster[i_db]);
// 	}
// 	printf("]\n");
	
	//	revisar cc
	int test1 = test_cc(cluster,*cluster_length);
	if(test1){
		//	for other nodes (recursividad) ( Intersection N(V Node_cluster) \ {cluster} == Intersection N(V Node_cluster))
		int candidates[Array_Size]={0};
		int candidates_len = 0;
		int i=0, test2 =0;
		int node = 0;
		
		//		Intersections
		node = cluster[0];
		candidates_len = adj_list_sizes[node];
		copy_array(adj_list[node],candidates,candidates_len);
		for (i=1;i < *cluster_length;i++){
			node =cluster[i];
			array_intersect_self(candidates,&candidates_len,adj_list[node],adj_list_sizes[node]);
		}
		
		//		recursividad
// 		printf("n candidates = %d \n",candidates_len);
		for(i=0;i< candidates_len;i++){
			node =candidates[i];
			test2 += conected_component(cluster,cluster_length,node,score_th);
		}
		
		//	Imprime si es necesario
		if(test2 == 0){
			
			//	Comprueba que sea nueva
			int test_in_class = in_class(cluster, *cluster_length);
			if(test_in_class == 2){
				fprintf(Log,"Error in test class\n");
				return 1;
			}
			if(test_in_class==0){
				//	Si no esta la clase entonces imprime
				print_cluster(cluster, *cluster_length,score_th);
				if(Classes_size>=Classes_max_size){
					fprintf(Log,"Class overload\n");
					return 1; 
				}
				//	Añade la clase (in_class ya lo hace)
				Classes_size++;
			}
		}
		
		remove_in_order_array(cluster,cluster_length,node);
		return 1;
	}else{
		//	dead end of search
		remove_in_order_array(cluster,cluster_length,node);
		return 0;
	}
}

//	Set operations

int in_class(int *cluster,int cluster_length){
	
	if(Classes_size >= Array_Size){return 2;}
	
	//	Obtain line_ID for every element of cluster
	int cluster_ID[Clases_adj_max]={0};
	int i =0, j= 0;
	for (i=0;i<cluster_length;i++){
		// 		printf("%u:",cluster[i]);
		cluster_ID[i]=line_ID[cluster[i]];
		// 		printf("%u\t",cluster_ID[i]);
	}
	// 	printf("\n");
	
	//	Test if the cluster is a sub set of any class
	for(i=0;i< Classes_size;i++){
		
		int test_subeq = is_subeq(cluster_ID,cluster_length,Classes[i],Classes_sizes[i]);
		if(test_subeq){
			//	It is a sub set of the i-th class
			return 1;
		}
	}
	
	//	If not, then add
	copy_array(cluster_ID,Classes[Classes_size],cluster_length);
	Classes_sizes[Classes_size]=cluster_length;
	
	return 0;
	
}

int is_subeq(int *Arr1,int size1, int *Arr2, int size2){
	int i=0,j=0;
	if(size1 >= Array_Size || size2 >= Array_Size){
		fprintf(Log,"Err is subeq in sizes \n");
		return 1;
	}
	while(i< size1 && j<size2){
		if(Arr1[i]>Arr2[j]){
			j++;
		}else if(Arr1[i]==Arr2[j]){
			i++;
			j++;
		}else if(Arr1[i]<Arr2[j]){
			return 0;
		}
	}
	if(i==size1){
		return 1;
	}
	return 0;
	
}

//	Funciones básicas

int max(int a , int b){
	if (a>= b){return a;}
	return b;
}

int min(int a , int b){
	if (a < b){return a;}
	return b;
}

//	Listas de adjacencia

int add_edge_adj_list( int nodo_a, int nodo_b){
	
	if(nodo_a >= Array_Size || nodo_b >= Array_Size ){return 1; }
	int i=0;
	
	// añade arista a -> b
	add_to_order_array( adj_list[nodo_a], &adj_list_sizes[nodo_a] ,  nodo_b);
	
	// añade arista b -> a
	add_to_order_array( adj_list[nodo_b], &adj_list_sizes[nodo_b] ,  nodo_a);
	
	return 0;
}
