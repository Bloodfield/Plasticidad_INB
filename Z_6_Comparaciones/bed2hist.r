
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1){
	stop("Only an output name should be provided (file extencion will be ignored).n",call.=FALSE)
}
name=args[1];

#	Variables

img_width=5000;
img_height=350;

#	Read Files
bed_table = read.table("stdin");

centromere_coords=paste(path.expand("~"),"/bin/centromere_coords.bed",sep="")
centromere_table = read.table(centromere_coords);
# centromere_table = read.table("/home/bloodfield/W_Root/1_Scy/Biología_Computer_Bioinformática/Plasticidad/Z_6_Comparaciones/centromere_coords.bed");
row.names(centromere_table) <- centromere_table[,1];
centromere_table <- centromere_table[,2-3]

Chr_end = paste(path.expand("~"),"/bin/human.hg38.genome",sep="")
Chr_end_table = read.table(Chr_end);
# Chr_end_table = read.table("/home/bloodfield/W_Root/1_Scy/Biología_Computer_Bioinformática/Plasticidad/Z_6_Comparaciones/human.hg38.genome");
row.names(Chr_end_table) <- Chr_end_table[,1]

#	For every Chromosome
chr_list = bed_table[,1];
chr_list = unique(chr_list);

#	nombre y guardado de imagenes
image_name = sub('\\..*$', '', name);
image_name = paste(image_name,".pdf",sep="");
# pdf(image_name, width = img_width, height = img_height);
pdf(image_name);
for (Chromosome in chr_list){
	
	#	Make compact data
	temp_chr_table_x = bed_table[bed_table[,1]==Chromosome,2];
	temp_chr_table_y = bed_table[bed_table[,1]==Chromosome,3];
	chr_table_x = c();
	chr_table_y_max = c();
	chr_table_y_avr = c();
	chr_table_y_min = c();
	step = temp_chr_table_x[length(temp_chr_table_x)] / img_width;
	for (i in 1:(img_width-1)) {
		chr_table_x = c(chr_table_x,(i-0.5)*step);
		segment_interval = (temp_chr_table_x >= (i-1)*step & temp_chr_table_x < i*step);
		segment = temp_chr_table_y[segment_interval];
		segment_data_avr = 0
		segment_data_min = 0
		segment_data_max = 0
		if(length(segment)>0){
# 			print(temp_chr_table_y[segment_interval])
			segment_data_max = max(segment);
			
			#	Promedio real
			inicio = chr_table_x[1]-1;
			fin = chr_table_x[length(chr_table_x)];
			factor_avr =  length(chr_table_x) / (fin - inicio);
			segment_data_avr = mean(segment)*factor_avr;
			
			#	Minimo real
			segment_data_min = 0;
			if (length(chr_table_x) == (fin - inicio)){
				segment_data_min = min(segment);
				segment_data_avr = mean(segment);
			}
		}
		chr_table_y_avr = c(chr_table_y_avr,segment_data_avr);
		chr_table_y_max = c(chr_table_y_max,segment_data_max);
		chr_table_y_min = c(chr_table_y_min,segment_data_min);
	}
	chr_table_x = c(chr_table_x,(img_width-0.5)*step,Chr_end_table[Chromosome,2]);
	print((img_width-0.5)*step)
	segment_interval = (temp_chr_table_x >= (img_width-1)*step & temp_chr_table_x <= img_width*step);
	segment = temp_chr_table_y[segment_interval];
	segment_data_avr = 0
	segment_data_min = 0
	segment_data_max = 0
	if(length(segment)>0){
# 			print(temp_chr_table_y[segment_interval])
		segment_data_max = max(segment);
		
		#	Promedio real
		inicio = chr_table_x[1]-1;
		fin = chr_table_x[length(chr_table_x)];
		factor_avr =  length(chr_table_x) / (fin - inicio);
		segment_data_avr = mean(segment)*factor_avr;
		
		#	Minimo real
		segment_data_min = 0;
		if (length(chr_table_x) == (fin - inicio)){
			segment_data_min = min(segment);
			segment_data_avr = mean(segment);
		}
	}
	chr_table_y_avr = c(chr_table_y_avr,segment_data_avr,0);
	chr_table_y_max = c(chr_table_y_max,segment_data_max,0);
	chr_table_y_min = c(chr_table_y_min,segment_data_min,0);
	
	
	
	#	Get Chromosome divitions
	
	height = max(chr_table_y_max);
	line = c(0,height);
	
	#	Histograma
	plot(	chr_table_x,
		chr_table_y_max,
		type="h",
		col="orangered3",
		main=Chromosome, 
		xlab="Coordenada", 
		ylab="Frecuencia de deleciones recurrentes",
		xlim=c(0, Chr_end_table[Chromosome,2]), 
		ylim=c(0, height),
		lwd = 0.5,
		axes=TRUE);
	Step_text = paste("Step = ",round(step)," bases",sep="");
	mtext(Step_text, side=3)
	lines(	chr_table_x,
		chr_table_y_avr,
		type="h",
		col="mediumseagreen",
		lwd = 0.5);
	lines(	chr_table_x,
		chr_table_y_min,
		type="h",
		col="blue",
		lwd = 0.5);
# 	axis(side=1,pos=Chr_end_table[Chromosome,2])
# 	axis(side=2,pos=height)
	
	#	Separaciones
	#		Mínimo
	div_coord = 0;	# cambiar por las coordenadas
	division_X=c(div_coord,div_coord);
	lines(division_X,line,type="l",lty=2,col="red");
	
	#		Máximo
	div_coord = Chr_end_table[Chromosome,2];	# cambiar por las coordenadas
	division_X=c(div_coord,div_coord);
	lines(division_X,line,type="l",lty=2,col="red");
	
	#		Centrómeros
	# for div_coord in div_list[Chromosome]
	if (any(row.names(centromere_table) == Chromosome)){
		div_coord = centromere_table[Chromosome,1];	# cambiar por las coordenadas
		division_X=c(div_coord,div_coord);
		lines(division_X,line,type="l",lty=2,col="red");
		
		div_coord = centromere_table[Chromosome,2];	# cambiar por las coordenadas
		division_X=c(div_coord,div_coord);
		lines(division_X,line,type="l",lty=2,col="red");
	}
}
dev.off();
warnings()
