#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1){
	stop("Only an output name should be provided (file extencion will be ignored).n",call.=FALSE)
}
name=args[1];

data = read.table("stdin");

pdf(name);
h <- hist(t(data),breaks=59);
plot(h)
dev.off();
