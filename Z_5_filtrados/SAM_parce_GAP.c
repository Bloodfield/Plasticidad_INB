//	RMFilter in C

/*	TODO
 * 	Mejorar Caché
 * 	Hacer guía
 * 
 */
#include <stdio.h>
#include <string.h>

#define Log_name "Log.txt"
#define INPUT_SIZE 4000
int clear_char(char *string, int size);
int get_line(char *string,int *size);
int get_info(char *line, int size,char *Name, char *Chromosome, char *CIGAR, unsigned *inicio );
int call_GAP(FILE *output,char *CIGAR,unsigned coordenada, unsigned C_min, char *Chromosome, char *Label);
int contained(char element,char *C_List,int size);
	
int main(int argc, char *argv[]){
	
	FILE *Log;
	Log = fopen(Log_name,"a");
	
	if (argc != 3){
		fprintf(Log,"Se tiene que escribir:\n");
		fprintf(Log,"\t Nombre base sin espacios \n");
		fprintf(Log,"\t Número de bases mínima para considerar un flanco \n");
		fprintf(Log,"$ SAM_parce_GAP nombre n_flanco  \n");
		fprintf(Log,"Emeplo:\n");
		fprintf(Log,"\n");
		fprintf(Log,"$ SAM_parce_GAP SRR3290534 20\n");
		fprintf(Log,"\n");
		fprintf(Log,"\n");
		fprintf(Log,"Este archivo el nombre base es \"SRR3290534\"\n");
		fprintf(Log,"Con el cual se llamarán los archivos de salida al añadirl el cromosoma\n");
		fprintf(Log,"\n");
		fprintf(Log,"\n");
		fclose(Log);
		return 1;
	}
	
	//	Obtiene argumentos
	
	char Base_name[300]={0};
	unsigned int n_flanco = 10;
	fprintf(Log,"Name : %s \n",argv[1]);
	strcpy(Base_name,argv[1]);
	sscanf (argv[2],"%u",&n_flanco);
	
	//	PARCEO de GAPS
	
	char input[INPUT_SIZE]={0};
	int size = 0;
	while (get_line(input,&size)!= EOF){
// 		printf("%s\n",input);
		
		if (input[0]!='@'){
			
			//	Get SAM info
			char Label[200]={0};
			char Chromosome[200]={0};
			char CIGAR[200]={0};
			unsigned inicio=0;
			int test = get_info(input, size,Label, Chromosome, CIGAR, &inicio );
			
			//	Open output file
			char out_name[300]={0};
			FILE *output;
			strcpy(out_name, Base_name);
			strcat(out_name,"_");
			strcat(out_name,Chromosome);
			strcat(out_name,"_");
			strcat(out_name,"splitreads.bed");
			
			output = fopen(out_name,"a");
			
			//	Parceo de CIGAR
			call_GAP(output,CIGAR,inicio,n_flanco, Chromosome, Label);
			
			fclose(output);
			
			clear_char(input,4000);
			clear_char(Label,200);
			clear_char(Chromosome,200);
			clear_char(CIGAR,200);
		}
	}
	

	fclose(Log);
	return 0;
}


//	Obtiene una linea de stdin
int get_line(char *string,int *size){
	*size = 0;
	char c;
	c = getchar();
	while (c  != '\n' && c != EOF){
		if (*size >= INPUT_SIZE){
			fprintf(Log,"Error : get_line : 0 : To long line %s\n",string);
			
		}
		string[*size]=c;
		(*size)++;
		c = getchar();
	}
	
	return c;
}

//	Maquina de estados para la lectura de la línea
int get_info(char *line, int size,char *Name, char *Chromosome, char *CIGAR, unsigned *inicio ){
	char state = 'A';
	char c;
	unsigned input_i = 0;
	int i=0;
	for (i=0;i < size; i++){
		c=line[i];
		switch(state){
			case 'A':
				if(c== '@'){
					state = 'B';
				}else{
					state = 'C';
					Name[input_i]=c;
					input_i++;
				}
				break;
			case 'B':
				if(c== '\n'){
					(*inicio) = 0;
					state = 'A';
				}
				break;
			case 'C':
				if(c!='\t'){
					Name[input_i]=c;
					input_i++;
				}else{
					input_i=0;
					state = 'D';
				}
				break;
			case 'D':
				if(c=='\t'){
					state = 'E';
				}
				break;
			case 'E':
				if(c!='\t'){
					Chromosome[input_i]=c;
					input_i++;
				}else{
					input_i=0;
					state = 'F';
				}
				break;
			case 'F':
				if(c!='\t'){
					(*inicio)=((*inicio)*10)+c-'0';
				}else{
					state = 'G';
				}
				break;
			case 'G':
				if(c=='\t'){
					state = 'H';
				}
				break;
				
			case 'H':
				if(c!='\t'){
					CIGAR[input_i]=c;
					input_i++;
				}else{
					input_i=0;
					state = 'B';
				}
				break;
				
		}
	}
	return 0;
}

//	Maquina de estados para el CIGAR
int call_GAP(FILE *output,char *CIGAR,unsigned coordenada, unsigned C_min,  char *Chromosome, char *Label){
	
	//	Variables principales
	
	char C_List[]={'I','S','M','=','X'};
	char G_List[]={'D','N'};
	char O_List[]={'H','P'};
	
	unsigned	C1	= 0;
	unsigned	DC2	= 0;
	unsigned	C2	= 0;
	unsigned	G	= 0;
	unsigned	inicio	= 0;
	unsigned	fin	= 0;
	unsigned	i	= 0;
	unsigned temp = 0;
	
	
	//	Completa C2 en base a los match que hay en el CIGAR
	
	while (CIGAR[i]!= 0 && i < 200){
		
		char c = CIGAR[i];
		
		if( c >= '0' && c <= '9'){
			temp = temp*10 +c -'0';
		}else if (contained(c, C_List,5)){
			C2 +=  temp;
			temp = 0;
			
		}else{
			temp = 0;
		}
		i++;
	}
	
	//	Maquina de estados
	
	
	char	state	= 'A';
	i=0;
	temp = 0;
	char c;
	while(CIGAR[i]!=0 && i < 200){
		//	Stand By
		c = CIGAR[i];
		if( c >= '0' && c <= '9'){
			temp = temp*10 +c -'0';
		}else{
			switch (state){
			case 'A':
				if(contained(c, C_List,5)){
					DC2 +=  temp;
					temp = 0;
					state = 'B';
				}else{
					temp=0;
				}
				break;
		
		//	Clevarage 1
			case 'B':
				if(contained(c,C_List,5)){
					DC2 +=  temp;
					temp = 0;
				}else if(contained(c,G_List,3)){
					G +=  temp;
					state = 'C';
					temp = 0;
					C2 -= DC2;
					C1 += DC2;
				}
				break;
		
		//	Gap
			case 'C':
				if(contained(c,G_List,3)){
					G += temp;
					temp = 0;
				}else if(contained(c,C_List,5)){
					DC2 +=   temp;
					temp = 0;
					state = 'D';
				}
				break;
		
		//	Clevarage 2	// puede dividirce entere B y C
			case 'D':
				if(contained(c,C_List,5)){
					DC2 +=  temp;
					temp = 0;
				}else if(contained(c,G_List,3)){
					
					//	Test -> print
					
					if(C1 >= C_min && C2 >= C_min && G > 0){
						inicio	= coordenada + C1;
						fin	= inicio + G;
						fprintf(output, "%s\t%u\t%u\t%s\n",Chromosome,inicio,fin,Label);
					}
					
					//	RESET
					C1	+= DC2;
					C2	-= DC2;
					DC2	= 0;
					G	= temp;
					inicio	= 0;
					fin	= 0;
					state	= 'C';
					temp = 0;
				}
				break;
		
		//	FSM Err
			default: fprintf(Log,"Error in CIGAR %s\n",CIGAR);
			}
		}
		i++;
	}
	
	//	Ultimo paso
	//	Test -> print
	if(C1 >= C_min && C2 >= C_min && G > 0){
		inicio	= coordenada + C1;
		fin	= inicio + G;
		fprintf(output, "%s\t%u\t%u\t%s\n",Chromosome,inicio,fin,Label);
	}
	
	return 0;
}

//	Revisa se el caracter "element" se encuentra en C_List
int contained(char element,char *C_List,int size){
	int i = 0;
	for (i =0; i<size; i++){
		if (element == C_List[i])
			return 1;
	}
	return 0;
}

//	Limpia string
int clear_char(char *string, int size){
	int i=0;
	for(i=0; i< size; i++){
		string[i]=0;
	}
	return 0;
}
