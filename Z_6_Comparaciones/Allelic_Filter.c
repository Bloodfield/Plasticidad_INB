//	RMFilter in C

/*	TODO
 * 	Hacer recvisión de flancos derecho e izquierdo, que en efecto sea izquierdo y derecho
 * 	Hacer guía
 * 	Funcionalizal main
 */
#include <stdio.h>

//	0_Core
#define Log_name "Log.txt"

int print_help(FILE *Log);

//	Z_Work

int main(int argc, char *argv[]){
	
	//	Head
	
	FILE *Log;
	Log = fopen(Log_name,"a");
	char *deleciones;
	char *flanco_izquierdo;
	char *flanco_derecho;
	float th=0.0;
		
	if (argc != 5){
		print_help(Log);
		fclose(Log);
		return 1;
	}
	
	fprintf(Log,"Name : %s \n",argv[1]);
	fprintf(Log,"Izquierdo : %s \n",argv[2]);
	fprintf(Log,"Derecho : %s \n",argv[3]);
	deleciones=argv[1];
	flanco_izquierdo=argv[2];
	flanco_derecho=argv[3];
	sscanf(argv[4],"%f",&th);
	
	//	Abrir y probar archivos
	
	FILE *deleciones_fh;
	FILE *flanco_izquierdo_fh;
	FILE *flanco_derecho_fh;
	deleciones_fh = fopen(deleciones,"r");
	flanco_izquierdo_fh = fopen(flanco_izquierdo,"r");
	flanco_derecho_fh = fopen(flanco_derecho,"r");
	
	int ret =0;
	if(!deleciones_fh ){
		fprintf(Log,"Archivo : %s no existe\n",deleciones);
		ret= 1;
	}
	if(!flanco_izquierdo_fh ){
		fprintf(Log,"Archivo : %s no existe\n",flanco_izquierdo);
		ret= 1;
	}
	if(!flanco_derecho_fh ){
		fprintf(Log,"Archivo : %s no existe\n",flanco_derecho);
		ret= 1;
	}
	if (ret){
		fclose(Log);
		return 1;
	}
	//	Algoritmo BEGIN
	while(!feof(deleciones_fh)){
		unsigned int inicio=0, fin=0;
		unsigned int cover_score_left=0, cover_score_right=0, cover_score_del=0;
		unsigned int lenght=0;
		float landscape_coverage=0.0, del_coverage=0.0;
		float ratio=0.0;
		char Chr[500]={0};
		
		if(fscanf(flanco_izquierdo_fh,"%u\n",&cover_score_left)!= 1){
			fprintf(Log,"Archivo : %s is not a BED format\n",flanco_izquierdo);
			fclose(Log);
			return 1;
		}
		
		if(fscanf(flanco_derecho_fh,"%u\n",&cover_score_right)!= 1){
			fprintf(Log,"Archivo : %s is not a BED format\n",flanco_derecho);
			fclose(Log);
			return 1;
		}
		if(fscanf(deleciones_fh,"%s\t%u\t%u\t%u\n",Chr,&inicio,&fin,&cover_score_del)!= 4){
			fprintf(Log,"Archivo : %s is not a BED format\n",deleciones);
			fclose(Log);
			return 1;
		}
		
		lenght = fin -inicio;
// 		printf("lenght = %u\n",lenght);
		landscape_coverage=(cover_score_left+cover_score_right)/40.0;
// 		printf("land cov = %0.6f\n",landscape_coverage);
		del_coverage=cover_score_del/(1.0*lenght);
// 		printf("del cov = %0.6f\n",del_coverage);
		ratio=del_coverage/landscape_coverage;
// 		printf("ratio = %0.6f\n",ratio);
		
		if((ratio < 1-th && ratio >0.5+th) || (ratio <0.5-th && ratio >th)){
			printf("%s\t%u\t%u\t%.5f\n",Chr,inicio,fin,ratio);
		}
	}
	
	//	Algoritmo END
	fclose(deleciones_fh);
	fclose(flanco_izquierdo_fh);
	fclose(flanco_derecho_fh);
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
