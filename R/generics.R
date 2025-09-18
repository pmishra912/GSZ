setGeneric("mGSZm",function(object, rnaseq=FALSE,comparisons=NULL,other.methods=FALSE,cor.method="pearson",min.cl.sz=5,max.cl.sz=500,pre.var=0,wgt1=0.2,wgt2=0.5,var.constant=10,perm.number=200) standardGeneric("mGSZm"))




setGeneric("exprdat", function(object) standardGeneric("exprdat"))
setGeneric("genesets", function(object) standardGeneric("genesets"))
setGeneric("sampleLabels", function(object) standardGeneric("sampleLabels"))


setGeneric("exprdat<-", function(object,value) standardGeneric("exprdat<-"))
setGeneric("genesets<-", function(object,value) standardGeneric("genesets<-"))
setGeneric("sampleLabels<-", function(object,value) standardGeneric("sampleLabels<-"))


### Generics for sub-classes mGSZResultAll, mGSZResultAll

setGeneric("topResults", function(object, direction = c("mixed","up","down"),order.by="P-value", n=10, method=NULL, ...) standardGeneric("topResults"))
