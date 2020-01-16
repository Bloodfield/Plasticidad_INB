#	HEAD

#	Run:
#	cp ${origin}/Graficador_KaryploteR.r .
#	cp ${origin}/Intersect_Recurrent_Cancer.bed .
#	Rscript Graficador_KaryploteR.r

#	Estructura:
#		HEAD
#		INIT
#		FUNCTIONS
#		MAIN

#	INIT
library("GenomicRanges")
library("karyoploteR")

# Constantes
chromosomes <- c("chr1",
	"chr2",
	"chr3",
	"chr4",
	"chr5",
	"chr6",
	"chr7",
	"chr8",
	"chr9",
	"chr10",
	"chr11",
	"chr12",
	"chr13",
	"chr14",
	"chr15",
	"chr16",
	"chr17",
	"chr18",
	"chr19",
	"chr20",
	"chr21",
	"chr22",
	"chrX",
	"chrY")
pdf_label <- c("[1 b, 10 b)",
	"[10 b, 100 b)",
	"[100 b, 1 Kb)",
	"[1 Kb, 10 Kb)",
	"[10 Kb, 100 Kb)",
	"[100 Kb, 1 Mb)",
	"[1 Mb, 10 Mb)",
	"[10 Mb, 100 Mb)",
	"[100 Mb, 500 Mb)")
pdf_name <- c("deletions_diez",
	"deletions_cien",
	"deletions_k",
	"deletions_diez_k",
	"deletions_cien_k",
	"deletions_m",
	"deletions_diez_m",
	"deletions_cien_m",
	"deletions_quinientos_m" )

#	Cargar datos
Deletions_table <- read.table("Deletions_Uniq_Recurrent.bed")
cancer <- read.table("cancer_higligth_chr.txt")
cancer_GR <- GRanges(seqnames = cancer$V1, ranges = IRanges(start = cancer$V2,end = cancer$V3))

#	Requerimientos y parámetros del plot
ymax <- 20
color <- "#853232"
color_cancer <- "#409478"
pp <-  getDefaultPlotParams(plot.type=2)
pp$data1height <- 100
pp$data2height <- 100
pp$topmargin <- 20
pp$bottommargin <- 10
pp$rightmargin <- 0.05
pp$data1outmargin <- 1
pp$data2outmargin <- 1
pp$ideogramheight <- 10
pp$data1inmargin <- 5
pp$data2inmargin <- 10

#	FUNCTIONS

#	Impresión
Plot_karyotyper <- function(temp,interval,name){
	
	chr <- unique(as.character(seqnames(temp)))

	for(i in chr){
		if(i %in% chromosomes){
			name_temp <- paste(name,"_",as.character(i),".pdf",sep="")
			pdf(name_temp,width = 12 , height = 3)
			
			kp <- plotKaryotype(plot.type=2,chromosomes=i, plot.params = pp)
			
			kpRect(kp,data = cancer_GR , y0=0, y1=1, col=color_cancer,border=NA)
			kpPlotCoverage(kp, data=temp,ymax=ymax, col=color)
			kpAxis(kp, ymax=ymax, cex=0.8)
			kpAbline(kp, h=4, lty=2, ymax=ymax)
			
			kpRect(kp,data = cancer_GR , y0=0, y1=1,col=color_cancer,border=NA, data.panel=2)
			kpPlotRegions(kp, data=temp, data.panel=2, col=color, layer.margin = 0.1, num.layers=20)
# 			kpAxis(kp, ymax=ymax, cex=0.8, data.panel=2)
			
			kpAddMainTitle(kp, main=interval)
			dev.off()
		}
	}
}


#	MAIN
#	Graficar para cada cromosoma en cada orden de magnitud

for(i in 1:9){
	interval_regs <- (Deletions_table$V3-Deletions_table$V2 >= 10^(i-1)) & (Deletions_table$V3-Deletions_table$V2 < 10^i)
	Del_temp = Deletions_table[ interval_regs , ]
	deletions <- GRanges(
		seqnames = Del_temp$V1,
		ranges = IRanges(start = Del_temp$V2,end = Del_temp$V3) )
	Plot_karyotyper(c(deletions),pdf_label[i],pdf_name[i])
}
