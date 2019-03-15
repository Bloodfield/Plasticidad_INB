
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1){
	stop("Only an output name should be provided (file extencion will be ignored).n",call.=FALSE)
}
name=args[1];

#	Variables

img_width=700;
img_height=350;

#	Read Files
bed_table = read.table("stdin");

centromere_coords=paste(path.expand("~"),"/bin/centromere_coords.bed",sep="")
centromere_table = read.table(centromere_coords);
row.names(centromere_table) <- centromere_table[,1];
centromere_table <- centromere_table[,2-3]

Chr_end = paste(path.expand("~"),"/bin/human.hg38.genome",sep="")
Chr_end_table = read.table(Chr_end);
row.names(Chr_end_table) <- Chr_end_table[,1]
# Chr_end_table <- Chr_end_table[,2]

#	For every Chromosome
chr_list = bed_table[,1];
chr_list = unique(chr_list);

#	nombre y guardado de imagenes
image_name = sub('\\..*$', '', name);
image_name = paste(image_name,".pdf",sep="");
pdf(image_name, width = img_width, height = img_height);

for (Chromosome in chr_list){
	
	#	Get Histogram values
	chr_table_x = bed_table[bed_table[,1]==Chromosome,2];
	chr_table_y = bed_table[bed_table[,1]==Chromosome,3];
	
	#	Get Chromosome divitions
	
	height = max(chr_table_y);
	line = c(0,height);
	
	#	Histograma
	plot(	chr_table_x,
		chr_table_y,
		type="h",
		col="blue",
		main=Chromosome, 
		xlab="Coordenada", 
		ylab="Cobertura",
		xlim=c(0, Chr_end_table[Chromosome,2]), 
		ylim=c(0, height));
# 		axes=FALSE);
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
