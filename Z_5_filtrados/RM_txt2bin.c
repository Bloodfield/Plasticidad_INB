//	RMFilter bin reference creation

//	Recibe el chr.txt > chr.bin
//	Reporta proceso en stdout

#include <stdio.h>
#include <string.h>

#define Log_name "Log.txt"

int makebin(FILE *txt_in, FILE *bin_out, FILE *Log);

int main (int argc, char *argv[]){
	
	
	char *rmsk_list;
	
	FILE *Log;
	Log = fopen(Log_name,"a");
	
	if (argc == 2){
		printf("Base name list is : %s \n",argv[1]);
		rmsk_list=argv[1];
	}else{
		printf("One reference list must be supplied in a file \n");
		printf("$ RMake rmsk_List.txt  \n");
		printf("\n");
		printf("\n");
		printf("\n");
		printf("Donde la lista tiene los nombres sin extención.\n Ejemplo:  \n");
		printf("$ head -5 rmsk_List.txt  \n");
		printf("$ rmsk_chr10_GL383545v1_alt  \n");
		printf("$ rmsk_chr10_GL383546v1_alt  \n");
		printf("$ rmsk_chr10_KI270824v1_alt  \n");
		printf("$ rmsk_chr10_KI270825v1_alt  \n");
		printf("$ rmsk_chr10  \n");
		return 1;
	}
	
	FILE *List;
	char in_name[500] = "";
	char out_name[500] = "";
	
	List = fopen(rmsk_list,"r");
	
	//	Por cada uno de los nombres en la lista, hacer conversión
	while (fscanf(List, "%s", in_name) == 1){
		
		//	Nombres de archivos de entrada y salida
		FILE *txt_in, *bin_out;
		strcpy(out_name,in_name);
		strcat(in_name,".txt");
		strcat(out_name,".bin");
		
		printf("%s\t>\t%s\n",in_name,out_name);
		
		//	Abrir
		txt_in = fopen(in_name,"r");
		bin_out= fopen(out_name,"w");
		
		//	Crea binario
		makebin( txt_in, bin_out,  Log);
		
		//	Cerrar
		fclose(txt_in);
		fclose(bin_out);
	}
	fclose(Log);
	fclose(List);
	return 0;
}


int makebin(FILE *txt_in, FILE *bin_out, FILE *Log){
	char c ;
	int tab_count = 0;
	long unsigned int coord_count = 0;
	char state = '0';
	//	Lee todo el archivo de entrada, y obtiene los valores numéricos de la posición 6 y 7 por cada línea
	while ((c = fgetc(txt_in)) != EOF){
		switch(state){
			//	wait tab
			case '0':
				if (c== '\t'){
					tab_count++;
					if (tab_count>=6){
						tab_count = 0;
						state = '1';
						
					}
				}
				break;
			//	Write number a
			case '1':
				if (c <= '9' && c >= '0' ){
					coord_count= coord_count*10+c-'0';
				}else if(c=='\t'){
					state = '2';
					fwrite(&coord_count, sizeof(long unsigned int), 1, bin_out);
					coord_count = 0;
				}
				else{
					fprintf(Log, "Error in c");
					state = 'z';
				}
				break;
			//	Write number b
			case '2':
				if (c <= '9' && c >= '0' ){
					coord_count= coord_count*10+c-'0';
				}else if(c=='\t'){
					state = 'z';
					fwrite(&coord_count, sizeof(long unsigned int), 1, bin_out);
					coord_count = 0;
				}
				else{
					fprintf(Log, "Error in c");
					state = 'z';
				}
				break;
			//	wait end of line
			case 'z':
				if (c=='\n'){state='0';}
			
		}
	}
	return 0;
}
