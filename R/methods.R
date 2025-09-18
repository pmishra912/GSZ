
setMethod("mGSZm",
          "GSZ",
          function(object,rnaseq=FALSE,comparisons=NULL,other.methods=FALSE,cor.method="pearson",min.cl.sz=5,max.cl.sz=500,pre.var=0,wgt1=0.2,wgt2=0.5,var.constant=10,perm.number=200){
          res <- .mGSZm(expr.data=exprdat(object),gene.sets=genesets(object),sample.labels=sampleLabels(object),rnaseq,comparisons,other.methods,cor.method="pearson",min.cl.sz,max.cl.sz,pre.var,wgt1,wgt2,var.constant,perm.number)
          return(res)
          }
)


setMethod("exprdat", "GSZ", function(object) object@exprdat)

setMethod("genesets", "GSZ", function(object) object@genesets)

setMethod("sampleLabels", "GSZ", function(object) object@sampleLabels)

setValidity("GSZ",function(object){
    msg <- NULL
    valid <- TRUE
    if(is.list(genesets(object))){
        gsets <- .toMatrix(exprdat(object),genesets(object))
        
        if(nrow(exprdat(object)) != nrow(gsets)){
            valid <- FALSE
            msg <- c(msg, "Number of genes in expression data and gene set data must be identical.")
        }
    }
    else
    if(nrow(exprdat(object)) != nrow(genesets(object))){
        valid <- FALSE
        msg <- c(msg, "Number of genes in expression data and gene set data must be identical.")
    }
})



setMethod("exprdat<-",signature(object="GSZ",value="matrixORdata.frame"),
function(object, value) {
    object@exprdat <- value
    if (validObject(object))
    return(object)
})


setMethod("genesets<-",signature(object="GSZ",value="matrixORdata.frameORlist"),
function(object, value) {
    object@genesets <- value
    if (validObject(object))
    return(object)
})

setMethod("sampleLabels<-",signature(object="GSZ",value="numericORcharacter"),
function(object, value) {
    object@sampleLabels <- value
    if (validObject(object))
    return(object)
})


#### Methods for sub-class mGSZResult and mGSZResultAll ###


setMethod("topResults","mGSZResult",function(object,direction=c("mixed","up","down"),order.by="P-value",n=10){
    res <- object@mGSZ
    
    if(order.by=="none"){
       return(res)
    }
    
   direction <- match.arg(direction)
    
    if(direction=="mixed"){
        res <- res
    }

    if(direction=="up"){
        res <- res[which(res[,6]=="up"),]
    }
    
    if(direction=="down"){
        res <- res[which(res[,6]=="down"),]
    }
    
    ind <- grep(order.by,colnames(res))
    
    if(ind==7|ind==8){
        res <- res[order(res[,ind]),]}
    else{
        res <- res[order(res[,ind],decreasing=T),]
    }
    return(res[1:n,])
    })




###############################################################################

setMethod("topResults","mGSZResultAll",function(object,direction=c("mixed","up","down"),order.by="P-value",n=10,method=NULL){
    direction <- match.arg(direction)
    if(is.null(method)){
        res <- object@mGSZ
    }
    else{
        res <- slot(object,method)
    }
    
    if(order.by=="none"){
        return(res)
    }
    
    
    if(direction=="mixed"){
        res <- res
    }

    if(direction=="up"){
        res <- res[which(res[,6]=="up"),]
    }
    
    if(direction=="down"){
        res <- res[which(res[,6]=="down"),]
    }
    
    
        ind <- grep(order.by,colnames(res))
        if(ind==7|ind==8){
            res <- res[order(res[,ind]),]}
        else{
            res <- res[order(res[,ind],decreasing=T),]
        }
    return(res[1:n,])
})


#### Methods for sub-class mGSZmResult###

setMethod("topResults","mGSZmResult", function(object,direction=c("mixed","up","down"),order.by="P-value",n=10,comparison,method=NULL){
    direction <- match.arg(direction)
    all_res <- object@results
    comp <- object@AllComparisons # "a-b" "a-c"
    
    
    
    if(length(comparison)==1){
    if(is.null(method)){
        res2 <- all_res[[comparison]]$mGSZm}
    else{
        res1 <- all_res[[comparison]]
        res2 <- res1[[method]]}
    
    if(order.by=="none"){
        return(res2)
    }
    
    if(direction=="mixed"){
        res <- res2}
    if(direction=="up"){
        res <- res2[which(res2[,6]=="up"),]}
    if(direction=="down"){
        res <- res2[which(res2[,6]=="down"),]}
    ind <- grep(order.by,colnames(res))
    if(ind==7|ind==8){
        res <- res[order(res[,ind]),]}
    else{
        res <- res[order(res[,ind],decreasing=T),]}
        
        if(n==Inf){
            return(res)
        }
        else{
            return(res[1:n,])
        }
    }
    
    
    else{### from here for multiple results
    
    noComps <- length(comparison)
    resultsSep <- list()
    for(i in 1:noComps){
        
        if(is.null(method)){
            res2 <- all_res[[comparison[i]]]$mGSZm}
        else{
            res1 <- all_res[[comparison[i]]]
            res2 <- res1[[method]]}
        if(direction=="mixed"){
            res <- res2}
        if(direction=="up"){
            res <- res2[which(res2[,6]=="up"),]}
        if(direction=="down"){
            res <- res2[which(res2[,6]=="down"),]}
        
        if(order.by=="none"){
            res <- res
        }
        
        else{
        ind <- grep(order.by,colnames(res))
        if(ind==7|ind==8){
            res <- res[order(res[,ind]),]
        }
        
        else{
           res <- res[order(res[,ind],decreasing=T),]
        }}
        
        
        
        if(n==Inf){
        
        # resultsSep[[i]] <- cbind(Contrast = comparison[i], res[,c(1,4,6,7)])
        resultsSep[[i]] <- cbind(Contrast = comparison[i], res)
        }
        
        else{
            
            #resultsSep[[i]] <- cbind(Contrast = comparison[i], res[1:n,c(1,4,6,7)])
            resultsSep[[i]] <- cbind(Contrast = comparison[i], res[1:n,])
        }
    }
        resultsSep <- do.call("cbind",resultsSep)
        rownames(resultsSep) <- NULL
        return(resultsSep)
    }

})


