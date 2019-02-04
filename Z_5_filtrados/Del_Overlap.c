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
#define Array_Size 1200
#define Classes_max_size 850
#define Clases_adj_max 1000
#define Adj_list_max_size 850
#define Str_len 300

unsigned de_queue(unsigned *array,unsigned end);
unsigned queue(unsigned *array,unsigned element, unsigned end, FILE *Log);
int percentage_overlap(unsigned reference , unsigned query ,unsigned *flanco_A,unsigned *flanco_B);
int err_message(FILE *Log);
int Del_Overlap(FILE *In_file,FILE *Log,unsigned overlap_th, unsigned score_th,char *Chr);
int read_line(FILE *In_file, unsigned *flanco_A, unsigned *flanco_B,unsigned *line_ID, unsigned *score, unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned *fin,unsigned overlap_th, FILE *Log);
unsigned in_class(unsigned *cluster,unsigned cluster_length, unsigned *line_ID, unsigned Classes[][Clases_adj_max],unsigned *Classes_sizes, unsigned Classes_size);
unsigned max(unsigned a , unsigned b);
int print_cluster(unsigned *flanco_A, unsigned *flanco_B, unsigned *factor, unsigned fin, unsigned *cluster, unsigned cluster_length, char *Chr);
int conected_component(unsigned *cluster,unsigned *cluster_length,unsigned node, unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes);
int copy_array(unsigned *array_a,unsigned *array_b,unsigned size);
int array_and_self(unsigned *array_a,unsigned *array_b,unsigned size);
int contract_to(unsigned *neighbor,unsigned fin,unsigned *cluster,unsigned *cluster_length);
unsigned add_element(unsigned *flanco_A, unsigned *flanco_B,unsigned *line_ID, unsigned line_number, unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned *fin,unsigned overlap_th,unsigned FA, unsigned FB, FILE *Log);
unsigned complete_data(FILE *In_file, unsigned *flanco_A, unsigned *flanco_B,unsigned *line_ID, unsigned *factor,  unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned *fin,unsigned overlap_th, FILE *Log);
unsigned remove_line_matrix(unsigned matrix[][Array_Size], unsigned n, unsigned lim_line, unsigned lim_row);
unsigned remove_row_matrix(unsigned matrix[][Array_Size], unsigned n, unsigned lim_line, unsigned lim_row);
unsigned recorrer_array(unsigned *array , unsigned size, unsigned n);
unsigned empty_array(unsigned *array, unsigned size);
unsigned remove_node_adj_list(unsigned adj_list[][Adj_list_max_size], unsigned *adj_list_sizes,unsigned node, unsigned fin);
unsigned add_to_order_array(unsigned *array, unsigned *len , unsigned elem);
unsigned array_intersect_self(unsigned *self_array,unsigned *array_len, unsigned *second,unsigned second_len);
unsigned add_edge_adj_list(unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned nodo_a, unsigned nodo_b);
unsigned is_subeq(unsigned *Arr1,unsigned size1, unsigned *Arr2, unsigned size2);
unsigned clear_unsigned_array(unsigned *Array, unsigned fin);
	
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
	
	// Cierre de archivos
	fclose(In_file);
	fclose(Log);
	
	if (Overflow==1){return 1;}
	
	
	return 0;
}

/*
 * 	PROGRAMAS DE SISTEMA
 * 	
 * 	DeQueue
 * 	Queue
 * 	percentage_overlap
 * 	err_message
 * 	Del_Overlap
 * 	empty_array
 * 	remove_line_matrix
 * 	remove_row_matrix
 * 	recorrer_array
 * 	complete_data
 * 	read_line
 * 	add_element
 * 	in_class
 * 	max
 * 	print_cluster
 * 	conected_component
 * 	copy_array
 * 	array_and_self
 * 	contract_to
 */


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
// 	printf("ref \t %u\t%u\nque\t%u\t%u\n",FAR,FBR,FAQ,FBQ);
	if (FBQ < FAR || FAQ > FBR){return 0;}
	unsigned divisor = FBR - FAR;
	unsigned temp= 0;
	if (FAQ < FAR){
		if(FBQ < FBR){
			temp = (FBQ-FAR )*100;
// 			printf("\ttemp1=%u\n",temp);
		}else{
			return 100;
		}
	}else{
		if(FBQ < FBR){
			temp = (FBQ-FAQ )*100;
// 			printf("\ttemp3=%u\n",temp);
		}else{
			temp = (FBR-FAQ )*100;
// 			printf("\ttemp4=%u\n",temp);
		}
	}
// 	printf("temp=%u \tdivisor = %u\n",temp,divisor);
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
	fprintf(Log,"Ejemplo:\n");
	fprintf(Log,"\n");
	fprintf(Log,"$ Del_Overlap $HOME/Deletions_Chr1.bed 50 5 Chr1\n");
	fprintf(Log,"\n");
	fprintf(Log,"El output se reporta en standard output \n");
	fprintf(Log,"\n");
	return 0;
}

int Del_Overlap(FILE *In_file,FILE *Log,unsigned overlap_th, unsigned score_th,char *Chr){
	
	unsigned Classes[Classes_max_size][Clases_adj_max]={0};
	unsigned Classes_sizes[Clases_adj_max];
	
	unsigned Classes_size=0;
	unsigned next_line=0;
	unsigned Overflow = 0;
	
	while(Overflow < 1){
		
// 		printf("Echo 2 \n");
// 		printf("Elems %u\n",fin);
		unsigned flanco_A[Array_Size]={0};
		unsigned flanco_B[Array_Size]={0};
		unsigned adj_list[Array_Size][Adj_list_max_size]={0};
		unsigned adj_list_sizes[Array_Size]={0};
		unsigned factor[Array_Size]={0};
		unsigned line_ID[Array_Size]={0};
		
		unsigned fin=0;
		
		
		//	Coloca la linea de lectura correspondiente
		fseek(In_file,next_line,SEEK_SET);
		
		//	Completa las lecturas con los datos del archivo si el archivo no ha fallado o terminado
		if(Overflow ==0){
			Overflow=complete_data(In_file, flanco_A, flanco_B,line_ID,factor,adj_list, adj_list_sizes, &fin, overlap_th,Log);
// 			printf("Echo 2.1 \n");
		}
		//	Si algo falló, regresa el fallo
		if(Overflow==2){
			return 1;
		}
		
		unsigned i = 0, j=0;
		
// 		for(i=0; i< fin;i++){
// 			printf("FA=%u\tFB=%u\tID=%u\tfactor=%u\n",flanco_A[i],flanco_B[i],line_ID[i],factor[i]);
// 		}
		
		//	Obten el conjunto de componentes conexas
		unsigned cluster[Array_Size]={0};
		unsigned cluster_length=0;
		conected_component(cluster,&cluster_length,0,adj_list, adj_list_sizes);
		
		unsigned test_in_class = in_class(cluster, cluster_length, line_ID,  Classes,Classes_sizes, Classes_size);
		
		if(test_in_class == 2){return 1;}
		if(test_in_class==0){
			//	Si no esta la clase entonces imprime
			print_cluster(flanco_A, flanco_B, factor, fin, cluster, cluster_length,Chr);
			if(Classes_size>=Classes_max_size){
				fprintf(Log,"Class overload\n");
// 				printf("Class overload\n");
				return 1; 
			}
			//	Añade la clase
			Classes_size++;
			}
		
		//	Purgar Clases
		
		//		Quitar linea actual de Clases
// 		printf("Purge next_line = %u\n",next_line);
		for(i=0;i<Classes_size;i++){
			for (j=0;j<Classes_sizes[i] && Classes[i][j] < next_line; j++){}
// 			printf("Purge j1 = %u \n",j);
			if(Classes[i][j] == next_line && Classes_sizes[i]>0 ){
				unsigned lim=Classes_sizes[i]-1;
				for (j;j<lim;j++){
					Classes[i][j] = Classes[i][j+1];
				}
				Classes[i][j] = 0;
				Classes_sizes[i]--;
			}
		}
		//		Eliminar clases vacías
		for(i=0;i<Classes_size && Classes_sizes[i]!=0;i++){}	//	Busca una linea vacía
		if(Classes_sizes[i]==0 && i<Classes_size){
			unsigned offset=1;
			//	Realiza corrimiento de datos
			while(Classes_sizes[i+offset]==0 && i+offset<Classes_size){
				offset++;
			}
			for(i;i<Classes_size-offset;i++){
				Classes_sizes[i]=Classes_sizes[i+offset];
				for(j=0;j<Classes_sizes[i];j++){
					Classes[i][j]=Classes[i+offset][j];
				}
				while(Classes_sizes[i+offset]==0 && i+offset< Classes_size-1){
					offset++;
				}
				if(Classes_sizes[i+offset]==0 && i+offset==Classes_size-1){
					offset++;
				}
			}
			
			//	Vuelve 0 las últimas
			for(i;i<Classes_size;i++){
				for (j;j<Classes_sizes[i];j++){
					Classes[i][j]=0;
				}
				Classes_sizes[i]=0;
			}
			
			//	Actualiza Classes_size
			Classes_size -= offset;
		}
		
		//	Point to next line in file
		if(fin > 1){
			next_line = line_ID[1];
		}else if(Overflow==0){
			next_line = ftell(In_file);
		}
// 		printf("next_line=%u\n",next_line);
		
		
		
		if (Overflow == 2){return 1;}
		
	}
		
	return 0;
	
}

unsigned clear_unsigned_array(unsigned *Array, unsigned fin){
	int i = 0;
	for (i=0;i<fin;i++){
		Array[i]=0;
	}
	return 0;
}

unsigned remove_node_adj_list(unsigned adj_list[][Adj_list_max_size], unsigned *adj_list_sizes, unsigned node, unsigned fin){
	
	unsigned i=0, j=0, next=0, fin_1 = fin -1;
	
	//	Recorre lista de adyacencia desde el nodo indicado
	if (fin >= Array_Size){return 1;}
	for (i=node; i< fin_1; i++){
		next = i+1;
		unsigned lim= adj_list_sizes[next];
		if (lim >= Array_Size){return 1;}
		for(j=0;j<lim;j++){
			adj_list[i][j]=adj_list[next][j];
		}
		adj_list_sizes[i]=adj_list_sizes[next];
	}
	for(j=0;j<adj_list_sizes[fin_1];j++){
		adj_list[i][j]=0;
	}
	adj_list_sizes[fin_1]=0;
	
	// reducir índices de los nodos afectados
	for (i=0; i< fin_1; i++){
		unsigned lim= adj_list_sizes[i];
		for (j=0;j<lim && adj_list[i][j] < node; j++){}
		if (adj_list[i][j]==node && lim >0){
			// En caso de encontrar el número eliminado, también quitarlo y recorrer
			lim--;
			for (j;j<lim;j++){
				adj_list[i][j] = adj_list[i][j+1]-1;
			}
			adj_list[i][j]=0;
			adj_list_sizes[i]--;
		}else{
			for (j;j<lim;j++){
				adj_list[i][j] --;
			}
		}
	}
	return 0;
	
}
unsigned empty_array(unsigned *array, unsigned size){
	unsigned i=0;
	if(size >= Array_Size){return 2;}
	for (i=0;i<size;i++){
		if(array[i]!=0){return 0;}
	}
	return 1;
}

unsigned remove_line_matrix(unsigned matrix[][Array_Size], unsigned n, unsigned lim_line, unsigned lim_row){
	
	if (n >= Array_Size || n >= lim_line){return 1;}
	unsigned lim = lim_line-1;
	unsigned i=0;
	
	for (i=n; i< lim; i++){
		copy_array(matrix[i+1],matrix[i],lim_row);
	}
	for (i=0;i< lim_row;i++){
		matrix[lim][i]=0;
	}
	return 0;
}

unsigned remove_row_matrix(unsigned matrix[][Array_Size], unsigned n, unsigned lim_line, unsigned lim_row){
	
	if (n >= Array_Size || n >= lim_line){return 1;}
	unsigned lim = lim_row-1;
	unsigned i=0;
	
	for (i=n; i< lim; i++){
		recorrer_array(matrix[i],lim_row,n);
	}
	return 0;
}

unsigned recorrer_array(unsigned *array, unsigned size, unsigned n){
	
	if(n >= size || n>= Array_Size ){return 1;}
	unsigned i=0;
	unsigned lim = size-1;
	for(i=n;i< lim;i++){
		array[i]=array[i+1];
	}
	array[lim]=0;
	return 0;
}

unsigned complete_data(FILE *In_file, unsigned *flanco_A, unsigned *flanco_B,unsigned *line_ID, unsigned *factor,  unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned *fin,unsigned overlap_th, FILE *Log){
	
	unsigned Overflow=read_line(In_file, flanco_A, flanco_B,line_ID,factor, adj_list,adj_list_sizes, fin, overlap_th,Log);
	
	unsigned FB0= flanco_B[0];
	unsigned FA0= flanco_A[0];
	//	Condition 1
	unsigned lim_buffer = FA0 +((FB0-FA0)*overlap_th/100)+1;	// +1 es por el redondeo
	unsigned FAX= 0;
	while(lim_buffer>FAX && Overflow == 0){
// 			printf("Echo 3 \n");
		//	Obten la información de la linea
		Overflow=read_line(In_file, flanco_A, flanco_B,line_ID,factor, adj_list,adj_list_sizes, fin, overlap_th,Log);
		if (Overflow >1){
			fprintf(Log,"Overlflow in read_lin\n");
			return 2;
		}
		//	Revisa el orden
		if((*fin)>1){
		if(flanco_A[(*fin)-1] < flanco_A[(*fin)-2]){
			fprintf(Log,"Bed file not sorted\n");
// 				printf("Bedfile not sorted\n");
			err_message(Log);
			fclose(Log);
			return 2;
		}
		}
		FAX= flanco_A[(*fin)-1];
	}
// 		printf("Echo 3.f \n");
	return Overflow;
}

int read_line(FILE *In_file, unsigned *flanco_A, unsigned *flanco_B,unsigned *line_ID, unsigned *factor, unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned *fin,unsigned overlap_th, FILE *Log){
	
	char Dummy1[Str_len]={0};
	char Dummy2[Str_len]={0};
	unsigned FAN=0,FBN=0;
	unsigned Overflow= 0;
	//	Leer linea
	unsigned line_number = ftell(In_file);
	if(fscanf(In_file,"%s\t%u\t%u\t%[^\n]s\n",Dummy1,&FAN,&FBN,Dummy2)!= 4){
		
		return 1;
	}
// 	printf("%s\t%u\t%u\t%s\n",Dummy1,FAN,FBN,Dummy2);
	
	
	if(*fin == 0){
		//	Lee primera linea
		Overflow = add_element(flanco_A, flanco_B, line_ID, line_number, adj_list, adj_list_sizes, fin,overlap_th,FAN, FBN, Log);
		factor[(*fin)-1]++;
// 		printf("Echo 4.1\tFAN=%u\tFBN=%u\tF_AN=%u\tF_BN=%u\n",FAN,FBN,flanco_A[0],flanco_B[0]);
	}else{
		//	Condiciones de addición:
		unsigned condition1=0,condition2=0;
		unsigned fin_1=(*fin)-1;
		long int FAP = flanco_A[fin_1];
		long int FBP = flanco_B[fin_1];
		long int FANI=FAN, FBNI=FBN;
		condition1 =((FBNI-FANI) > (((FBP-FAP)*overlap_th)/100) ); //	Evita añadir cosas muy pequeñas
		condition2 =(FBNI < ((((FBP*100)+(FANI*(overlap_th-100)))/(overlap_th))+1));	// evita añadir cosas muy grandes
		if(condition1 && condition2){
			
			if(FAP !=FAN || FBP !=FBN){
				Overflow = add_element(flanco_A, flanco_B,line_ID, line_number, adj_list, adj_list_sizes, fin,overlap_th,FAN, FBN, Log);
				fin_1=(*fin)-1;
	// 			printf("Echo 4.3\n");
			}
			factor[fin_1]++;
		}
// 		printf("Echo 4.2\n");
	}
	
	if (Overflow > 0){return  2;}
	
	return 0;
}

unsigned add_element(unsigned *flanco_A, unsigned *flanco_B,unsigned *line_ID,unsigned line_number,  unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned *fin,unsigned overlap_th,unsigned FA, unsigned FB, FILE *Log){
	if( queue(flanco_A,FA,*fin,Log)){return 1;}
	if(queue(flanco_B,FB,*fin,Log)){return 1;}
	if(queue(line_ID,line_number,*fin,Log)){return 1;}
	unsigned i = 0;
	for (i=0;i< *fin ; i++){
		int per_ida = percentage_overlap(i,*fin,flanco_A,flanco_B);
		int per_regreso = percentage_overlap(*fin,i,flanco_A,flanco_B);
// 		printf("ida = %d \t regreso = %d\n",per_ida,per_regreso);
		if(per_ida > overlap_th && per_regreso > overlap_th){
			if(adj_list_sizes[i] >= Adj_list_max_size||adj_list_sizes[*fin] >= Adj_list_max_size){
				fprintf(Log,"AjList overflow\n");
// 				printf("AjList overflow\n");
				return 1;
				
			}
			add_edge_adj_list( adj_list,adj_list_sizes, i,*fin);
// 			matrix[i][*fin]=matrix[*fin][i]=1;
		}
	}
	(*fin)++;
	return 0;
}

unsigned add_edge_adj_list(unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes, unsigned nodo_a, unsigned nodo_b){
	
	if(nodo_a >= Array_Size || nodo_b >= Array_Size ){return 1; }
	unsigned i=0;
	
	// añade arista a -> b
	add_to_order_array( adj_list[nodo_a], &adj_list_sizes[nodo_a] ,  nodo_b);
	
	// añade arista b -> a
	add_to_order_array( adj_list[nodo_b], &adj_list_sizes[nodo_b] ,  nodo_a);
	
	return 0;
}

unsigned add_to_order_array(unsigned *array, unsigned *len , unsigned elem){
	int i =0;
	for (i=0; i < *len && array[i]<elem;i++){}
	unsigned temp = elem;
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


//	Checks if the index lxuster "cluster" is in "Classes". If not, it's added
unsigned in_class(unsigned *cluster,unsigned cluster_length,unsigned *line_ID , unsigned Classes[][Clases_adj_max],unsigned *Classes_sizes, unsigned Classes_size){
	
	if(Classes_size >= Array_Size){return 2;}
	
	//	Obtain line_ID for every element of cluster
	unsigned cluster_ID[Clases_adj_max]={0};
	unsigned i =0, j= 0;
	for (i=0;i<cluster_length;i++){
// 		printf("%u:",cluster[i]);
		cluster_ID[i]=line_ID[cluster[i]];
// 		printf("%u\t",cluster_ID[i]);
	}
// 	printf("\n");
	
	//	Test if the cluster is a sub set of any class
	for(i=0;i< Classes_size;i++){
		
		unsigned test_subeq = is_subeq(cluster_ID,cluster_length,Classes[i],Classes_sizes[i]);
		if(test_subeq){
			//	It is a sub set of the i-th class
			return 1;
		}
	}
	
	//	If not, then add
	copy_array(cluster_ID,Classes[Classes_size],cluster_length);
	Classes_sizes[Classes_size]=cluster_length;
// // 	printf("classes size = %u\n",Classes_size);
// 	for (i=0;i<cluster_length;i++){
// // 		printf("%u:",cluster[i]);
// // 		printf("%u:",cluster_ID[i]);
// // 		Classes[Classes_size][i]=cluster_ID[i];
// // 		printf("%u\t",Classes[Classes_size][i]);
// 	}
// 	printf("\n");
	
	return 0;
	
}

unsigned is_subeq(unsigned *Arr1,unsigned size1, unsigned *Arr2, unsigned size2){
	unsigned i=0,j=0;
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

unsigned max(unsigned a , unsigned b){
	if (a>= b){return a;}
	return b;
}

int print_cluster(unsigned *flanco_A, unsigned *flanco_B, unsigned *factor, unsigned fin, unsigned *cluster, unsigned cluster_length, char *Chr){
	unsigned max_coord=flanco_A[0];
	unsigned i = 0, score = 0;
// 	printf("cluster _len= %u\n",cluster_length);
	for (i=0;i<cluster_length;i++){
		unsigned index = cluster[i];
		if (index >= fin){return -1;}
// 		printf("score= %u\n",factor[index]);
		score += factor[index];
		max_coord = max(max_coord,flanco_B[index]);
	}
	printf("%s\t%u\t%u\t%u\n",Chr,flanco_A[0],max_coord,score);
	
	return 0;
}

int conected_component(unsigned *cluster,unsigned *cluster_length,unsigned node, unsigned adj_list[][Adj_list_max_size],unsigned *adj_list_sizes){
	unsigned Overflow=0;
	
	unsigned neighbor[Adj_list_max_size]={0};
	unsigned ref[Adj_list_max_size]={0};
	unsigned ref_len=adj_list_sizes[node],neighbor_len=adj_list_sizes[node];
	copy_array(adj_list[0],neighbor,neighbor_len);
	copy_array(adj_list[0],ref,ref_len);
	
	add_to_order_array( neighbor, &neighbor_len , node);
// 	printf("Echo 5.1 \n");
	unsigned i = 0;
	for (i=0; i<= ref_len ; i++){
		unsigned index = ref[i];
		unsigned temp[Adj_list_max_size]={0};
		unsigned temp_len=adj_list_sizes[index];
		Overflow= copy_array(adj_list[index],temp,temp_len);
// 		printf("echo 5.2 %u\n",Overflow);
		Overflow+=add_to_order_array( temp, &temp_len , index);
// 		printf("echo 5.3 %u\n",Overflow);
		Overflow+=array_intersect_self(neighbor,&neighbor_len, temp, temp_len);
// 		printf("echo 5.4 %u\n",Overflow);
		unsigned id=0,jd=0;
// 		for (id = 0; id < temp_len;id++){
// 			printf("%u\t",temp[id]);
// 		}
// 		printf("\n");
// 		for (id = 0; id < neighbor_len;id++){
// 			printf("%u\t",neighbor[id]);
// 		}
// 		printf("\n");
	}
// 	printf("Echo 5.f \n");
	(*cluster_length)=neighbor_len;
	copy_array(neighbor,cluster,*cluster_length);
	
	return 0;
}

unsigned array_intersect_self(unsigned *self_array,unsigned *array_len, unsigned *second,unsigned second_len){
	
	if(*array_len>Adj_list_max_size || second_len >Adj_list_max_size){return 1;}
	unsigned temp[Adj_list_max_size]={0};
	unsigned temp_len=0;
	unsigned i=0,j=0;
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

int copy_array(unsigned *array_a,unsigned *array_b,unsigned size){
	unsigned i=0;
	for(i=0;i<size; i++){
		array_b[i]=array_a[i];
	}
// 	printf("Echo 6 \n");
	return 0;
}

int array_and_self(unsigned *array_a,unsigned *array_b,unsigned size){
	unsigned i=0;
	for(i=0;i<size; i++){
		array_a[i]= array_b[i] && array_a[i];
	}
	return 0;
}

int contract_to(unsigned *neighbor,unsigned fin,unsigned *cluster,unsigned *cluster_length){
	unsigned i=0;
	if(fin >= Array_Size){return -1;}
	(*cluster_length)=0;
	for(i=0;i<fin; i++){
		if(neighbor[i]){
			cluster[*cluster_length]=i;
			(*cluster_length)++;
		}
	}
	return 0;
}
