library(qtl)

##### Load genotype and phenotype file
geno= read.cross("csv",  file="Geno_BC2S3_Pheno_div461.csv",na.strings = "-",genotypes = c("AA", "AB","BB"), alleles = c("A", "B"), estimate.map = FALSE, BC.gen=2, F.gen=3)
##### jittermap function required 
geno <- jittermap(geno)
summary(geno)
##################################### calculate QTL conditional probability ############################################
geno <- calc.genoprob(geno, step=0.5, off.end=0, error.prob=0.0001, map.function="kosambi", stepwidth="fixed")

##################################### QTL mapper, including scanone and mqm ############################################
mapper= function(geno, pheno.col, threshold ) {

#################################################### primary single-QTL scanning #############################

out<- scanone (geno,pheno.col=pheno.col,method="hk")
ylimit= max(as.data.frame(out[,-(1:2)]))
ylimit


##### scanone results summary
single.summary=summary(out,threshold=threshold,lodcolumn=1) 
single.summary

######### output scanone results
write.csv(out, file=paste(pheno.col,"_scanone.csv",sep=""), quote=FALSE,row.names=T)
######### output scanone summarized results
write.csv(single.summary, file=paste(pheno.col,"_scanone_summary.csv",sep=""), quote=FALSE,row.names=T)
#################################################### end of single-QTL mapping #########################



#################################################### beginning of mqm mapping ##########################
if(dim(single.summary)[1] !=0 ){
looper=1
##### make initial QTL object for fitting multiple QTL model
chr= single.summary$chr
pos= single.summary$pos
qtl= makeqtl(geno, chr=chr, pos=pos, what="prob")
qtl

while(looper ==1){

##### construct formula of multiple QTL model
term= rownames(summary(qtl))
formula1= as.formula(paste("y ~ ", paste(term, collapse= "+")))
formula1

##### refine QTL model
qtl= refineqtl(geno, pheno.col=pheno.col, qtl=qtl, formula= formula1, method= "hk")
qtl

##### scan whole genome to check if there are additional significant QTL
out.aq <- addqtl(geno, pheno.col=pheno.col, qtl=qtl, formula=formula1, method="hk")
single.summary= summary(out.aq,threshold=threshold,lodcolumn=1 )
single.summary

if(dim(single.summary)[1] !=0 ){

##### if there does more significant QTLs are detected, add into previous QTL object
rqtl.added= addtoqtl(geno, qtl, single.summary$chr, single.summary$pos)
qtl= rqtl.added

looper=1

}else{
looper=0
}

}

##### plot final multiple QTL mapping results
pdf(paste(pheno.col,"_mqm.pdf",sep=""), width = 6, height=4)
plotLodProfile(qtl,main=pheno.col,ylab="LOD",cex=0.4,lwd=1)
dev.off()

##### output multiple QTL fitting results
out.mqm= fitqtl(geno, pheno.col=pheno.col, qtl=qtl, formula= formula1, method= "hk", dropone=T, get.ests=T)
summary(out.mqm)

write.table(summary(out.mqm)[[1]], file= paste(pheno.col,"_mqm_full.txt",sep=""), quote=FALSE, sep='\t', row.names=TRUE, col.names=T)
write.table(summary(out.mqm)[[2]], file= paste(pheno.col,"_mqm_dropone.txt",sep=""), quote=FALSE, sep='\t', row.names=TRUE, col.names=T)
write.table(summary(out.mqm)[[3]], file= paste(pheno.col,"_mqm_effect.txt",sep=""), quote=FALSE, sep='\t', row.names=TRUE, col.names=T)


######### Here is QTL scanning results with other QTL fixed as background
names(attributes(qtl))                                                                       
lodprof <- attr(qtl, "lodprofile")  #### results saved in the component lodprofile                                                                                    
lodprof[[1]] #### first QTL scanning results                                                                                                           
lodint(qtl, qtl.index=1)  ##### support interval for the first QTL                                                                                                                                                                                          

######## output scanning results for each QTL, saved in lodprof.txt
out.lodprof<-capture.output(lodprof)
cat(out.lodprof,file=paste(pheno.col,"_mqm_lodprof.txt",sep=""),sep="\n",append=TRUE)


######## output QTL support interval for each QTL, saved in qtlint.txt
out.int=NULL ###### The purpose here is to empty the possible previously existing file
cat(out.int,file=paste(pheno.col,"_mqm_qtlint.txt",sep=""),sep="\n")
for (i in 1:length(lodprof)){                                                                                                       
qtlint= lodint(qtl, qtl.index=i)
out.int<-capture.output(qtlint)
cat(out.int,file=paste(pheno.col,"_mqm_qtlint.txt",sep=""),sep="\n",append=TRUE)
}

}
}
######################################################### end of QTL mapper ######################



######################################################### analyze your traits ###############################

pheno.list=names(geno$pheno)
###### loop through your traits
for(i in 1:length(pheno.list)){
mapper(geno=geno, pheno.col=pheno.list[i], threshold=5)
}
######################################################### end of mapping #############################################
