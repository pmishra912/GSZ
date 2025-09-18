.mGSZm <- function(expr.data,rnaseq = FALSE,gene.sets,sample.labels,comparisons=NULL,other.methods=NULL,cor.method="pearson", min.cl.sz,max.cl.sz, pre.var,wgt1,wgt2,var.constant,perm.number){

    # rnaseq: TRUE for RNAseq data matrix
    
    start.val=5
    

   ########################################
   ## Testing data sizes ##
   ########################################
   
   list.of.packages <- c("ismev", "MASS","edgeR","DESeq2")
   new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
   if(length(new.packages)) install.packages(new.packages)



   suppressMessages(suppressWarnings(require(ismev)))
   suppressMessages(suppressWarnings(require(MASS)))
   suppressMessages(suppressWarnings(require(edgeR)))
   suppressMessages(suppressWarnings(require(DESeq2)))
   


    if(is.null(comparisons) & length(unique(sample.labels))==2){
     print("Gene set analysis of two-group gene expression data")
        result <- .mGSZ(expr.data=expr.data,rnaseq=rnaseq,gene.sets=gene.sets,sample.labels=sample.labels,min.cl.sz=min.cl.sz,max.cl.sz=max.cl.sz,other.methods=other.methods,pre.var=pre.var,wgt1=wgt1,wgt2=wgt2,var.constant=var.constant,perm.number=perm.number)
        
    }


  if(is.null(comparisons) & length(unique(sample.labels))!=2){
    print("Gene set analysis with continuous phenotype")
    result <- .mGSZ_cont(expr.data=expr.data,rnaseq=rnaseq,gene.sets=gene.sets,sample.labels=sample.labels,cor.method=cor.method,min.cl.sz=min.cl.sz,max.cl.sz=max.cl.sz,other.methods=other.methods,pre.var=pre.var,wgt1=wgt1,wgt2=wgt2,var.constant=var.constant,perm.number=perm.number)
    
 }


if(!is.null(comparisons) & length(unique(sample.labels))>2){
    print("Gene set analysis with > 2 groups. Advanced permutation technique will be used if the number of replicates is < 6.")
    result <- .mGSZmm(expr.data=expr.data,rnaseq=rnaseq,gene.sets=gene.sets,sample.labels=sample.labels,min.cl.sz=min.cl.sz,max.cl.sz=max.cl.sz,other.methods=other.methods,comparisons=comparisons, pre.var=pre.var,wgt1=wgt1,wgt2=wgt2,var.constant=var.constant,perm.number=perm.number)
    
}

return(result)
}


    ########################
    
    
    
.mGSZmm <- function(expr.data,rnaseq,gene.sets,sample.labels,comparisons,min.cl.sz,max.cl.sz,other.methods,pre.var,wgt1,wgt2,var.constant,perm.number){
    
    start.val=5
    
    # reordering expr.data and sample.labels
    ord_lbs <- order(sample.labels)
    sample.labels <- sample.labels[ord_lbs]
    expr.data <- expr.data[,ord_lbs]
    
    
    gene.sets <- .toMatrix(expr.data,gene.sets)
    expr.dat.sz <- length(unique(rownames(expr.data)))
    num.genes <- length(unique(rownames(gene.sets)))
    
    ##Removal of rows with no gene set members
    
    data <- .rm.rows.with.noSetMembers(expr.data,gene.sets)
    gene.sets <- data$gene.sets
    expr.data <- data$expr.data
    
    
    
    ##Removal of too small gene sets
    
    gene.sets <- .filter.genesets(gene.sets,min.cl.sz, max.cl.sz)
    geneset.classes <- colnames(gene.sets)
    num.classes <- ncol(gene.sets)
    geneset.class.sizes <- colSums(gene.sets)
    
    
    
    ##Calculation of differential gene expression scores
    
    tmp.expr.data = .des_perm4(expr.data, gene.sets, labels=sample.labels, rnaseq=rnaseq, comparisons=comparisons, perm.number)
    
    
    
    
    diff.expr.dat.sz <- dim(tmp.expr.data[[1]])
    
    
    ## mGSZ score for positive data ####
    
    pos.gsz <- list()
    pos.gsz.direction <- list()
    pos.gsa <- list()
    pos.allez <- list()
    upper.tail <- list()
    lower.tail <- list()
    
    if(other.methods){
        pos.ss <- list()
        pos.sum <- list()
        pos.wks <- list()
    }
    
    
    for(i in 1:length(comparisons)){
        res <- .mGSZ.test.score(expr.data=tmp.expr.data[[i]][,1], gene.sets, wgt1, wgt2, pre.var, var.constant, start.val)
        pos.gsz[[i]] <- res$gene.set.scores$mGSZ.scores
        pos.gsz.direction[[i]] <- res$gene.set.scores$mGSZ.direction
        
        upper.tail[[i]] <- res$gene.set.scores$upper.tail
        lower.tail[[i]] <- res$gene.set.scores$lower.tail
        
        
        pos.gsa[[i]] <- res$gene.set.scores$mGSA.scores
        pos.allez[[i]] <- res$gene.set.scores$mAllez.scores
        
        if(other.methods){
            pos.ss[[i]] <- .SS.test.score(expr.data=tmp.expr.data[[i]][,1], gene.sets)
            pos.sum[[i]] <- .sumTestscore(expr.data=tmp.expr.data[[i]][,1], gene.sets)
            pos.wks[[i]] <- .KS.score(expr.data=tmp.expr.data[[i]][,1], gene.sets)$KS.mod
        }
    }
    names(pos.gsz) <- comparisons
    names(pos.gsa) <- comparisons
    names(pos.allez) <- comparisons
    
    
    if(other.methods){
        names(pos.ss) <- comparisons
        names(pos.sum) <- comparisons
        names(pos.wks) <- comparisons
    }
    
    
    ## mGSZ for sample permutation data ###
    
    col.perm.gsz <- list()
    col.perm.gsa <- list()
    col.perm.allez <- list()
    
    if(other.methods){
        col.perm.ss <- list()
        col.perm.sum <- list()
        col.perm.wks <- list()
    }
    
    
    for(j in 1:length(comparisons)){
        
        col.perm.gsz[[j]] <- matrix(0, perm.number, num.classes)
        col.perm.gsa[[j]] <- matrix(0, perm.number, num.classes)
        col.perm.allez[[j]] <- matrix(0, perm.number, num.classes)
        
        if(other.methods){
            col.perm.ss[[j]] <- matrix(0, perm.number, num.classes)
            col.perm.sum[[j]] <- matrix(0, perm.number, num.classes)
            col.perm.wks[[j]] <- matrix(0, perm.number, num.classes)
        }
    }
    
    
    for(l in 1:length(comparisons)){
        for( k in 1:perm.number){
            res <- .mGSZ.test.score(tmp.expr.data[[l]][,1+k], gene.sets,wgt1, wgt2, pre.var, var.constant, start.val)
            col.perm.gsz[[l]][k,] <- as.numeric(res$gene.set.scores$mGSZ.score)
            col.perm.gsa[[l]][k,] <- as.numeric(res$gene.set.scores$mGSA.score)
            col.perm.allez[[l]][k,] <- as.numeric(res$gene.set.scores$mAllez.score)
            
            if(other.methods){
                col.perm.ss[[l]][k,] <- .SS.test.score(tmp.expr.data[[l]][,1+k], gene.sets)
                col.perm.sum[[l]][k,] <- .sumTestscore(tmp.expr.data[[l]][,1+k], gene.sets)
                col.perm.wks[[l]][k,] <- .KS.score(tmp.expr.data[[l]][,1+k], gene.sets)$KS.mod
            }
            
        }
    }
    
    ### p-value calculation for the gene set scores with sample permutation ###
    pvals.gsz <- list()
    pvals.gsa <- list()
    pvals.allez <- list()
    
    if(other.methods){
        pvals.ss <- list()
        pvals.sum <- list()
        pvals.wks <- list()
    }
    for(l in 1:length(comparisons)){
        pvals.gsz[[l]] <- .mGSZ.p.values(pos.gsz[[l]],col.perm.gsz[[l]])
        pvals.gsa[[l]] <- .mGSA.p.values(pos.gsa[[l]],col.perm.gsa[[l]])
        pvals.allez[[l]] <- .mAllez.p.values(pos.allez[[l]],col.perm.allez[[l]])
        
        
        if(other.methods){
            pvals.ss[[l]] <- .SS.p.values(pos.ss[[l]], col.perm.ss[[l]])
            pvals.sum[[l]] <- .mAllez.p.values(pos.sum[[l]], col.perm.sum[[l]])
            pvals.wks[[l]] <- .KS.p.values(pos.wks[[l]], col.perm.wks[[l]])
            
        }
        
    }
    
    names(pvals.gsz) <- names(pos.gsz)
    names(pvals.gsa) <- names(pos.gsa)
    names(pvals.allez) <- names(pos.allez)
    if(other.methods){
        names(pvals.ss) <- names(pos.ss)
        names(pvals.sum) <- names(pos.sum)
        names(pvals.wks) <- names(pos.wks)
    }
    
    
    
    
    if(rnaseq){
        vst <- varianceStabilizingTransformation(expr.data)
        expdat <- vst
    }
    
    else{
        expdat <- expr.data
    }
    
    #fc_scores <- .gene_sets_fc_score_mGSZm(expdat, comparisons, sample.labels, gene.sets)
    
    ##Output table for each of the methods
    
    result <- list()
    for(a in 1:length(comparisons)){
        #table.gsz <- data.frame(colnames(gene.sets)[pvals.gsz[[a]][[2]]], geneset.class.sizes[pvals.gsz[[a]][[2]]], upper.tail[[a]][pvals.gsz[[a]][[2]]], lower.tail[[a]][pvals.gsz[[a]][[2]]], pos.gsz[[a]][pvals.gsz[[a]][[2]]], pos.gsz.direction[[a]][pvals.gsz[[a]][[2]]], fc_scores[[a]][[2]][pvals.gsz[[a]][[2]]], fc_scores[[a]][[1]][pvals.gsz[[a]][[2]]], pvals.gsz[[a]][[1]],pvals.gsz[[a]][[3]],row.names=NULL)
        
        #colnames(table.gsz) <- c("Gene sets","No. of member genes","upper tail", "lower tail", "GSZ scores","direction","Mean FC","Mean absolute FC","P-value","FDR")
        
        table.gsz <- data.frame(colnames(gene.sets)[pvals.gsz[[a]][[2]]], geneset.class.sizes[pvals.gsz[[a]][[2]]], upper.tail[[a]][pvals.gsz[[a]][[2]]], lower.tail[[a]][pvals.gsz[[a]][[2]]], pos.gsz[[a]][pvals.gsz[[a]][[2]]], pos.gsz.direction[[a]][pvals.gsz[[a]][[2]]], pvals.gsz[[a]][[1]],pvals.gsz[[a]][[3]],row.names=NULL)
        
        colnames(table.gsz) <- c("Gene sets","No. of member genes","upper tail", "lower tail", "GSZ scores","direction","P-value","FDR")
        
        
        
        table.gsa <- data.frame(colnames(gene.sets)[pvals.gsa[[a]][[2]]], geneset.class.sizes[pvals.gsa[[a]][[2]]], pos.gsa[[a]][pvals.gsa[[a]][[2]]],  pvals.gsa[[a]][[1]],pvals.gsa[[a]][[3]],row.names=NULL)
        
        colnames(table.gsa) <- c("Gene sets","No. of member genes","mGSA scores","P-value","FDR")
        
        
        table.allez <- data.frame(colnames(gene.sets)[pvals.allez[[a]][[2]]],  geneset.class.sizes[pvals.allez[[a]][[2]]], pos.allez[[a]][pvals.allez[[a]][[2]]],pvals.allez[[a]][[1]],pvals.allez[[a]][[3]],row.names=NULL)
        
        colnames(table.allez) <- c("Gene sets","No. of member genes","mAllez scores","P-value","FDR")
        
        if(other.methods){
            
            table.ss <- data.frame(colnames(gene.sets)[pvals.ss[[a]][[2]]],  geneset.class.sizes[pvals.ss[[a]][[2]]], pos.ss[[a]][pvals.ss[[a]][[2]]],  pvals.ss[[a]][[1]],pvals.ss[[a]][[3]],row.names=NULL)
            colnames(table.ss) <- c("Gene sets","No. of member genes","SS scores","P-value","FDR")
            
            
            table.sum <- data.frame(colnames(gene.sets)[pvals.sum[[a]][[2]]],  geneset.class.sizes[pvals.sum[[a]][[2]]], pos.sum[[a]][pvals.sum[[a]][[2]]], pvals.sum[[a]][[1]],pvals.sum[[a]][[3]],row.names=NULL)
            colnames(table.sum) <- c("Gene sets","No. of member genes","SUM scores","P-value","FDR")
            
            table.wks <- data.frame(colnames(gene.sets)[pvals.wks[[a]][[2]]],  geneset.class.sizes[pvals.wks[[a]][[2]]], pos.wks[[a]][pvals.wks[[a]][[2]]],  pvals.wks[[a]][[1]],pvals.wks[[a]][[3]],row.names=NULL)
            colnames(table.wks) <- c("Gene sets","No. of member genes","wKS scores","P-value","FDR")
            
        }
        
        
        result[[a]] <- list(mGSZm = table.gsz, sample.labels = sample.labels, expr.data = expr.data, gene.sets = gene.sets,comparison=comparisons[a])
        
        if(other.methods){
            result[[a]] <- list(mGSZm = table.gsz, mGSAm = table.gsa, mALLEZm = table.allez, SSm = table.ss, SUMm = table.sum, wKSm = table.wks, sample.labels = sample.labels, expr.data = expr.data, gene.sets = gene.sets,comparison=comparisons[a])
        }
        
    }
    #names(result) <- paste("comparison",c(1:length(result)),sep="")
    names(result) <- comparisons
    out <- mGSZmResults(results=result,AllComparisons=comparisons)
    return(out)
}



############ mGSZ ###############

.mGSZ <- function(expr.data,rnaseq,gene.sets,sample.labels,min.cl.sz,max.cl.sz,other.methods,pre.var,wgt1,wgt2,var.constant,perm.number){

    start.val=5
    
    ord_lbs <- order(sample.labels)
    sample.labels <- sample.labels[ord_lbs]
    expr.data <- expr.data[,ord_lbs]
    
    
    
    
    
   ########################################
   ## Testing that data sizes are correct #
   ########################################
   expr.dat.sz <- dim(expr.data)
   gene.sets <- .toMatrix(expr.data,gene.sets)
   num.genes <- nrow(gene.sets)

   ##Remove rows with no gene set members
   data <- .rm.rows.with.noSetMembers(expr.data,gene.sets)
   expr.data <- data$expr.data
   gene.sets <- data$gene.sets
    
   ##Remove too small gene sets
   gene.sets <- .filter.genesets(gene.sets,min.cl.sz,max.cl.sz)
   geneset.classes <- colnames(gene.sets)
   num.classes <- dim(gene.sets)[2]
   num.genes <- dim(expr.data)[1]
   geneset.class.sizes <- colSums(gene.sets)
   
    ##Calculation of gene scores 
   tmp = .diffScore(expr.data, rnaseq, sample.labels, perm.number)
   tmp.expr.data = tmp$t_scores
   perm.number = tmp$perm.number
   diff.expr.dat.sz <- dim(tmp.expr.data)
	
	
	## mGSZ score for positive data ####
	pos.scores <- .mGSZ.test.score(expr.data=tmp.expr.data[,1], gene.sets, wgt1, wgt2, pre.var, var.constant, start.val)
    pos.mGSZ.scores <- pos.scores$gene.set.scores$mGSZ.scores
    pos.mGSZ.direction <- pos.scores$gene.set.scores$mGSZ.direction
    upper.tail <- pos.scores$gene.set.scores$upper.tail
    lower.tail <- pos.scores$gene.set.scores$lower.tail

  
    ## Gene set scores for positive data with other methods ##
    if(other.methods){
		pos.mAllez.scores <- pos.scores$gene.set.scores$mAllez.scores
		pos.mGSA.scores <- pos.scores$gene.set.scores$mGSA.scores
        pos.ss.scores <- .SS.test.score(tmp.expr.data[,1], gene.sets)
        pos.sum.scores <- .sumTestscore(tmp.expr.data[,1], gene.sets)
        pos.ks.scores <- .KS.score(tmp.expr.data[,1], gene.sets)$KS.org
        pos.wks.scores <- .KS.score(tmp.expr.data[,1], gene.sets)$KS.mod
        pos.ks.up <- .KS.score(tmp.expr.data[,1], gene.sets)$profile2

    }
	 
	## mGSZ for permutation data ###
	col.perm.mGSZ <- matrix(0, diff.expr.dat.sz[2]-1, num.classes)
	col.perm.mAllez <- col.perm.mGSZ
	col.perm.mGSA <- col.perm.mGSZ
    col.perm.WRS <- col.perm.mGSZ
    col.perm.SS <- col.perm.mGSZ
    col.perm.SUM <- col.perm.mGSZ
    col.ks <- col.perm.mGSZ
    col.wks <- col.perm.mGSZ
    

    for( k in 1:(diff.expr.dat.sz[2]-1)){
    	tmp <- .mGSZ.test.score(tmp.expr.data[,k+1], gene.sets,wgt1, wgt2, pre.var, var.constant, start.val)
    	col.perm.mGSZ[k,]<-tmp$gene.set.scores$mGSZ.scores

        if(other.methods){
    		col.perm.mAllez[k,] <- tmp$gene.set.scores$mAllez.scores
    		col.perm.mGSA[k,] <- tmp$gene.set.scores$mGSA.scores
            col.ks[k,] <- .KS.score(tmp.expr.data[,k+1], gene.sets)$KS.org
            col.wks[k,] <- .KS.score(tmp.expr.data[,k+1], gene.sets)$KS.mod
            col.perm.SS[k,] <- .SS.test.score(tmp.expr.data[,k+1], gene.sets)
            col.perm.SUM[k,] <- .sumTestscore(tmp.expr.data[,k+1], gene.sets)
        }
     }
     
     ### p-value calculation for the gene set scores with sample permutation ###

    mGSZ.p.vals.col.perm <- .mGSZ.p.values(pos.mGSZ.scores,col.perm.mGSZ)

    if(other.methods){
        ss.p.vals.col <- .SS.p.values(pos.ss.scores,col.perm.SS)
        sum.p.vals.col <- .mAllez.p.values(pos.sum.scores,col.perm.SUM)
        wKS.p.vals.col <- .KS.p.values(pos.wks.scores,col.wks)
        mGSA.p.vals.col.perm <- .mGSA.p.values(pos.mGSA.scores,col.perm.mGSA)
        mAllez.p.vals.col.perm <- .mAllez.p.values(pos.mAllez.scores,col.perm.mAllez)
    }


    
## calculation of mean FC and mean absolute FC for member genes of each gene sets


if(rnaseq){
    vst <- varianceStabilizingTransformation(expr.data)
    expdat <- vst
}


else{
    expdat <- expr.data
}

#fc_scores <- .expr.fc(expdat, sample.labels, gene.sets)
   
    
    
### preparing output table for each of the methods with column permutation ###

#mGSZ.table.col <- data.frame(colnames(gene.sets)[mGSZ.p.vals.col.perm$class.ind], geneset.class.sizes[mGSZ.p.vals.col.perm$class.ind],upper.tail[mGSZ.p.vals.col.perm$class.ind], lower.tail[mGSZ.p.vals.col.perm$class.ind], pos.mGSZ.scores[mGSZ.p.vals.col.perm$class.ind],pos.mGSZ.direction[mGSZ.p.vals.col.perm$class.ind], fc_scores$fc[mGSZ.p.vals.col.perm$class.ind], fc_scores$abs.fc[mGSZ.p.vals.col.perm$class.ind],mGSZ.p.vals.col.perm$EV.class,mGSZ.p.vals.col.perm$FDR,row.names=NULL)
#   colnames(mGSZ.table.col) <- c("Gene sets","No. of member genes","upper tail", "lower tail", "GSZ scores","direction","Mean FC","Mean absolute FC","P-value","FDR")

    mGSZ.table.col <- data.frame(colnames(gene.sets)[mGSZ.p.vals.col.perm$class.ind], geneset.class.sizes[mGSZ.p.vals.col.perm$class.ind],upper.tail[mGSZ.p.vals.col.perm$class.ind], lower.tail[mGSZ.p.vals.col.perm$class.ind], pos.mGSZ.scores[mGSZ.p.vals.col.perm$class.ind],pos.mGSZ.direction[mGSZ.p.vals.col.perm$class.ind],mGSZ.p.vals.col.perm$EV.class,mGSZ.p.vals.col.perm$FDR,row.names=NULL)
  colnames(mGSZ.table.col) <- c("Gene sets","No. of member genes","upper tail", "lower tail", "GSZ scores","direction","P-value","FDR")



    if(other.methods){
        mGSA.table.col <- data.frame(colnames(gene.sets)[mGSA.p.vals.col.perm$class.ind], geneset.class.sizes[mGSA.p.vals.col.perm$class.ind],pos.mGSA.scores[mGSA.p.vals.col.perm$class.ind],mGSA.p.vals.col.perm$EV.class,mGSA.p.vals.col.perm$FDR,row.names=NULL)
        colnames(mGSA.table.col) <- c("Gene sets","No. of member genes","mGSA scores","Mean FC","Mean absolute FC","P-value","FDR")
        
   
        mAllez.table.col <- data.frame(colnames(gene.sets)[mAllez.p.vals.col.perm$class.ind], geneset.class.sizes[mAllez.p.vals.col.perm$class.ind],pos.mAllez.scores[mAllez.p.vals.col.perm$class.ind],mAllez.p.vals.col.perm$NORM.class,mAllez.p.vals.col.perm$FDR,row.names=NULL)
        colnames(mAllez.table.col) <- c("Gene sets","No. of member genes","mAllez scores","Mean FC","Mean absolute FC","P-value","FDR")
    
        SUM.table.col <- data.frame(colnames(gene.sets)[sum.p.vals.col$class.ind],geneset.class.sizes[sum.p.vals.col$class.ind],pos.sum.scores[sum.p.vals.col$class.ind],sum.p.vals.col$NORM.class,sum.p.vals.col$FDR,row.names=NULL)
        colnames(SUM.table.col) <- c("Gene sets","No. of member genes","SUM scores","Mean FC","Mean absolute FC","P-value","FDR")
   
        SS.table.col <- data.frame(colnames(gene.sets)[ss.p.vals.col$class.ind], geneset.class.sizes[ss.p.vals.col$class.ind],pos.ss.scores[ss.p.vals.col$class.ind],ss.p.vals.col$EMP.class,ss.p.vals.col$FDR,row.names=NULL)
        colnames(SS.table.col) <- c("Gene sets","No. of member genes","SS scores","Mean FC","Mean absolute FC","P-value","FDR")
    
        New.KS.table.col <- data.frame(colnames(gene.sets)[wKS.p.vals.col$class.ind], geneset.class.sizes[wKS.p.vals.col$class.ind],pos.wks.scores[wKS.p.vals.col$class.ind],wKS.p.vals.col$EMP.class,wKS.p.vals.col$FDR,row.names=NULL)
        colnames(New.KS.table.col) <- c("Gene sets","No. of member genes","wKS scores","Mean FC","Mean absolute FC","P-value","FDR")
    }
    

  if(other.methods){
      out <- mGSZResultsAll(mGSZ=mGSZ.table.col,mGSA = mGSA.table.col, mAllez = mAllez.table.col,wKS = New.KS.table.col,SUM = SUM.table.col,SS = SS.table.col,sampleLabels=sample.labels,exprdat=expr.data,genesets=gene.sets)}
  
 else{
     out <- mGSZResults(mGSZ=mGSZ.table.col,sampleLabels=sample.labels,exprdat=expr.data,genesets=gene.sets)}
    return(out)
}

############
.diffScore <-
function(data, rnaseq, labels, perm.number) {

# data: expression data matrix
# labels: Vector of response values (example: 1,2)
# perms.number: Number of sample permutations

    
    dime2 <- dim(data)
  pit <- length(labels)
  
  if( !( is.array(perm.number) || is.matrix(perm.number)) ) {
    
  all_perms <- replicate(perm.number,sample(labels,pit,replace=FALSE))}
  unique.perm <- unique(t(all_perms))
  if(dim(unique.perm)[1]<perm.number){
    
        all_perms <- t(unique.perm)
  }
  dime  <- dim(all_perms)
  dime2 <- dim(data)
  diff_t = array(0, c(dime2[1], dime[2] + 1))
  diff_p = diff_t
  
### positive results ###

  length(labels)
  dim(data)
    
  if(rnaseq){
  print("RNAseq data inputed. This software expects raw counts!")
    dgelist <- DGEList(counts=data,group=labels) # creating edgeR's data object
    dgelist <- calcNormFactors(dgelist, method="TMM") # TMM normalization
    design <- model.matrix(~labels)
    y <- voom(counts=dgelist,design)
    fit1 <- lmFit(y, design)
        
  }
    
 else{
   mySet <- new("ExpressionSet",expr=as.matrix(data))
   design <- model.matrix(~labels)
   fit1 <- lmFit(mySet, design)
 }
    
  fit2<-eBayes(fit1)
  diff_t[,1] <- fit2$t[,2]
  diff_p[,1] <- fit2$p.value[,2]
  
## permutation results ###

  for(k in 1:dime[2]){
    design <- model.matrix(~all_perms[,k])  
    if(rnaseq){
        dgelist <- DGEList(counts=data,group=all_perms[,k]) # creating edgeR's data object
        dgelist <- calcNormFactors(dgelist, method="TMM") # TMM normalization
        y <- voom(counts=dgelist,design)
        fit1 <- lmFit(y, design)
    }
      
    else{
      fit1 <- lmFit(mySet, design)
    }
      
    fit2<-eBayes(fit1)
    diff_t[,k+1] <- fit2$t[,2]
    diff_p[,k+1] <- fit2$p.value[,2]
    }
  out <- list(p_values = diff_p, t_scores = diff_t, perms = unique.perm, perm.number=dim(unique.perm)[1])
  }

##################

.cbind_diff_dims <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

##########################################################################################################



.mGSZ_cont <- function(expr.data,rnaseq,gene.sets,sample.labels,cor.method,min.cl.sz,max.cl.sz,other.methods,pre.var,wgt1,wgt2,var.constant,perm.number){
    
    start.val=5
    
    
    ord_lbs <- order(sample.labels)
    sample.labels <- sample.labels[ord_lbs]
    expr.data <- expr.data[,ord_lbs]
    
    
    
    
    
    
    ########################################
    ## Testing that data sizes are correct #
    ########################################
    
    
    expr.dat.sz <- dim(expr.data)
    gene.sets <- .toMatrix(expr.data,gene.sets)
    num.genes <- nrow(gene.sets)
    
    ##Remove rows with no gene set members
    
    data <- .rm.rows.with.noSetMembers(expr.data,gene.sets)
    expr.data <- data$expr.data
    gene.sets <- data$gene.sets
    
    ##Remove too small gene sets
    
    gene.sets <- .filter.genesets(gene.sets,min.cl.sz,max.cl.sz)
    geneset.classes <- colnames(gene.sets)
    num.classes <- dim(gene.sets)[2]
    num.genes <- dim(expr.data)[1]
    geneset.class.sizes <- colSums(gene.sets)
    
    ##Calculation of gene scores
    
    tmp.expr.data <- .diffScoreCont(data=expr.data,pheno=sample.labels,cor.method="pearson",rnaseq=rnaseq,perm.number)
    
    
    diff.expr.dat.sz <- dim(tmp.expr.data)
    
    
    ## mGSZ score for positive data ####
    
    
    pos.scores <- .mGSZ.test.score(expr.data=tmp.expr.data[,1], gene.sets, wgt1, wgt2, pre.var, var.constant, start.val)
    
    pos.mGSZ.scores <- pos.scores$gene.set.scores$mGSZ.scores
    pos.mGSZ.direction <- pos.scores$gene.set.scores$mGSZ.direction
    
    upper.tail <- pos.scores$gene.set.scores$upper.tail
    lower.tail <- pos.scores$gene.set.scores$lower.tail
    
    
    ## Gene set scores for positive data with other methods ##
    
    if(other.methods){
        pos.mAllez.scores <- pos.scores$gene.set.scores$mAllez.scores
        pos.mGSA.scores <- pos.scores$gene.set.scores$mGSA.scores
        pos.ss.scores <- .SS.test.score(tmp.expr.data[,1], gene.sets)
        pos.sum.scores <- .sumTestscore(tmp.expr.data[,1], gene.sets)
        pos.ks.scores <- .KS.score(tmp.expr.data[,1], gene.sets)$KS.org
        pos.wks.scores <- .KS.score(tmp.expr.data[,1], gene.sets)$KS.mod
        pos.ks.up <- .KS.score(tmp.expr.data[,1], gene.sets)$profile2
        
    }
    
    
    ## mGSZ for permutation data ###
    
    
    col.perm.mGSZ <- matrix(0, diff.expr.dat.sz[2]-1, num.classes)
    col.perm.mAllez <- col.perm.mGSZ
    col.perm.mGSA <- col.perm.mGSZ
    col.perm.WRS <- col.perm.mGSZ
    col.perm.SS <- col.perm.mGSZ
    col.perm.SUM <- col.perm.mGSZ
    col.ks <- col.perm.mGSZ
    col.wks <- col.perm.mGSZ
    
    
    for( k in 1:(diff.expr.dat.sz[2]-1)){
        tmp <- .mGSZ.test.score(tmp.expr.data[,k+1], gene.sets,wgt1, wgt2, pre.var, var.constant, start.val)
        col.perm.mGSZ[k,]<-tmp$gene.set.scores$mGSZ.scores
        
        if(other.methods){
            col.perm.mAllez[k,] <- tmp$gene.set.scores$mAllez.scores
            col.perm.mGSA[k,] <- tmp$gene.set.scores$mGSA.scores
            col.ks[k,] <- .KS.score(tmp.expr.data[,k+1], gene.sets)$KS.org
            col.wks[k,] <- .KS.score(tmp.expr.data[,k+1], gene.sets)$KS.mod
            col.perm.SS[k,] <- .SS.test.score(tmp.expr.data[,k+1], gene.sets)
            col.perm.SUM[k,] <- .sumTestscore(tmp.expr.data[,k+1], gene.sets)
        }
    }
    
    
    
    
    ### p-value calculation for the gene set scores with sample permutation ###
    
    mGSZ.p.vals.col.perm <- .mGSZ.p.values(pos.mGSZ.scores,col.perm.mGSZ)
    
    if(other.methods){
        ss.p.vals.col <- .SS.p.values(pos.ss.scores,col.perm.SS)
        sum.p.vals.col <- .mAllez.p.values(pos.sum.scores,col.perm.SUM)
        wKS.p.vals.col <- .KS.p.values(pos.wks.scores,col.wks)
        mGSA.p.vals.col.perm <- .mGSA.p.values(pos.mGSA.scores,col.perm.mGSA)
        mAllez.p.vals.col.perm <- .mAllez.p.values(pos.mAllez.scores,col.perm.mAllez)
    }
    
    
    
    ## calculation of mean FC and mean absolute FC for member genes of each gene sets
    
    
    if(rnaseq){
        vst <- varianceStabilizingTransformation(expr.data)
        expdat <- vst
    }
    
    
    else{
        expdat <- expr.data
    }
    
    #fc_scores <- .expr.fc(expdat, sample.labels, gene.sets)
    
    
    
    ### preparing output table for each of the methods with column permutation ###
    
    mGSZ.table.col <- data.frame(colnames(gene.sets)[mGSZ.p.vals.col.perm$class.ind], geneset.class.sizes[mGSZ.p.vals.col.perm$class.ind],upper.tail[mGSZ.p.vals.col.perm$class.ind], lower.tail[mGSZ.p.vals.col.perm$class.ind],pos.mGSZ.scores[mGSZ.p.vals.col.perm$class.ind],pos.mGSZ.direction[mGSZ.p.vals.col.perm$class.ind],mGSZ.p.vals.col.perm$EV.class,mGSZ.p.vals.col.perm$FDR,row.names=NULL)
    colnames(mGSZ.table.col) <- c("Gene sets","No. of member genes","upper tail", "lower tail", "GSZ scores","direction","Mean FC","Mean absolute FC","P-value","FDR")
    
    if(other.methods){
        mGSA.table.col <- data.frame(colnames(gene.sets)[mGSA.p.vals.col.perm$class.ind], geneset.class.sizes[mGSA.p.vals.col.perm$class.ind],pos.mGSA.scores[mGSA.p.vals.col.perm$class.ind],mGSA.p.vals.col.perm$EV.class,mGSA.p.vals.col.perm$FDR,row.names=NULL)
        colnames(mGSA.table.col) <- c("Gene sets","No. of member genes","mGSA scores","Mean FC","Mean absolute FC","P-value","FDR")
        
        
        mAllez.table.col <- data.frame(colnames(gene.sets)[mAllez.p.vals.col.perm$class.ind], geneset.class.sizes[mAllez.p.vals.col.perm$class.ind],pos.mAllez.scores[mAllez.p.vals.col.perm$class.ind],mAllez.p.vals.col.perm$NORM.class,mAllez.p.vals.col.perm$FDR,row.names=NULL)
        colnames(mAllez.table.col) <- c("Gene sets","No. of member genes","mAllez scores","Mean FC","Mean absolute FC","P-value","FDR")
        
        SUM.table.col <- data.frame(colnames(gene.sets)[sum.p.vals.col$class.ind],geneset.class.sizes[sum.p.vals.col$class.ind],pos.sum.scores[sum.p.vals.col$class.ind],sum.p.vals.col$NORM.class,sum.p.vals.col$FDR,row.names=NULL)
        colnames(SUM.table.col) <- c("Gene sets","No. of member genes","SUM scores","Mean FC","Mean absolute FC","P-value","FDR")
        
        SS.table.col <- data.frame(colnames(gene.sets)[ss.p.vals.col$class.ind], geneset.class.sizes[ss.p.vals.col$class.ind],pos.ss.scores[ss.p.vals.col$class.ind],ss.p.vals.col$EMP.class,ss.p.vals.col$FDR,row.names=NULL)
        colnames(SS.table.col) <- c("Gene sets","No. of member genes","SS scores","Mean FC","Mean absolute FC","P-value","FDR")
        
        New.KS.table.col <- data.frame(colnames(gene.sets)[wKS.p.vals.col$class.ind], geneset.class.sizes[wKS.p.vals.col$class.ind],pos.wks.scores[wKS.p.vals.col$class.ind],wKS.p.vals.col$EMP.class,wKS.p.vals.col$FDR,row.names=NULL)
        colnames(New.KS.table.col) <- c("Gene sets","No. of member genes","wKS scores","Mean FC","Mean absolute FC","P-value","FDR")
    }
    
    
    if(other.methods){
        out <- mGSZResultsAll(mGSZ=mGSZ.table.col,mGSA = mGSA.table.col, mAllez = mAllez.table.col,wKS = New.KS.table.col,SUM = SUM.table.col,SS = SS.table.col,sampleLabels=sample.labels,exprdat=expr.data,genesets=gene.sets)}
    
    else{
        out <- mGSZResults(mGSZ=mGSZ.table.col,sampleLabels=sample.labels,exprdat=expr.data,genesets=gene.sets)}
    return(out)
}


############
.diffScoreCont <-
function(data, pheno, cor.method,rnaseq, perm.number) {
    
    # data: expression data matrix
    # pheno: Vector of continuous phenotype
    # perms.number: Number of sample permutations
    
    if(rnaseq){
        data <- varianceStabilizingTransformation(data)
    }
    
    perms <- replicate(perm.number,sample(c(1:ncol(data)),ncol(data),replace=FALSE))
    
    dime  <- dim(perms)
    dime2 <- dim(data)
    corr = array(0, c(dime2[1], dime[2] + 1))
    
    
    ### positive results ###
    corr[,1] <- apply(data,1,function(x) cor(x,pheno,method=cor.method))
    
    
    ## permutation results ###
    
    for(k in 1:dime[2]){
        corr[,k+1] <- apply(data[,perms[,k]],1,function(x) cor(x,pheno,method=cor.method))
        
    }
    return(corr)
}





#### sub-functions #####


### calculation of fold change of member genes ###
.gene_sets_fc_score_mGSZm <- function(expr.data, comparisons, sample.labels, gene.sets){
    fc <- list()
    for(j in 1:length(comparisons)){
        grps <- .getNumeric(comparisons[j])
        ind <- match(sample.labels,grps)
        ind.na <- which(is.na(ind))
        expr <- expr.data[,-ind.na]
        labels <- sample.labels[-ind.na]
        fc[[j]] <- .expr.fc(expr, labels, gene.sets)
    }
    return(fc)
}

## R code to calculate mean fold change/absolute fold change of expression of member genes
.expr.fc <- function(expr, labels, gene.sets){
    #expr.data: Subset of expr.data for a particular comparisons
    fc <- rep(NA,ncol(gene.sets))
    abs.fc <- fc
    grp1 <- table(labels)[[1]]
    grp2 <- table(labels)[[2]]
    for(i in 1:ncol(gene.sets)){
        ind <- which(gene.sets[,i]==1)
        exprdat <- expr[ind,]
        #fc[i] <- -mean(apply(exprdat[,1:grp1],1,mean)-apply(exprdat[,(grp1+1):(grp1+grp2)],1,mean))
        #abs.fc[i] <- mean(abs(apply(exprdat[,1:grp1],1,mean)-apply(exprdat[,(grp1+1):(grp1+grp2)],1,mean)))
        fc[i] <- mean(apply(exprdat[,(grp1+1):(grp1+grp2)],1,mean)-apply(exprdat[,1:grp1],1,mean)) # B-A form
        abs.fc[i] <- mean(abs(apply(exprdat[,(grp1+1):(grp1+grp2)],1,mean)-apply(exprdat[,1:grp1],1,mean)))
    }
    out <- list(abs.fc=abs.fc,fc=fc)
    return(out)
}

### R code to get numeric part from a string ###
.getNumeric <- function(string){
    str <- strsplit(string,"[^[:digit:]]")
    res <- as.numeric(unlist(str))
    sol <- res[!is.na(res)]
    return(sol)
}
#####################
###############
.toTable <- function(mGSZMobj,comparison,n){
    out <- mGSZMobj[[comparison]][[1]][1:n,]
    out <- out[order(out[,6],decreasing=F),]
    out
}

################

.rm.rows.with.noSetMembers <-
function(expr.data,gene.sets){
    tmp_sum <- apply(gene.sets, 1, sum)
    ind_wth_sets <- which(tmp_sum > 0)
    expr.data <- expr.data[ind_wth_sets,]
    gene.sets   <- gene.sets[ind_wth_sets,]
    out <- list(expr.data=expr.data, gene.sets=gene.sets)}

###############
## removes too small and too big gene sets ##
.filter.genesets <-
function(gene.sets,min.cl.sz, max.cl.sz){
    tmp_sum <- apply(gene.sets, 2, sum)
    ind_wth_sets <- which(tmp_sum > min.cl.sz & tmp_sum < max.cl.sz)
    gene.sets   <- gene.sets[,ind_wth_sets]
    out=gene.sets}

##################

.sumVarMean_calc <-
function(expr_data, gene.sets, pre.var){
    dim_sets <- dim(gene.sets)
    length_expr <- length(expr_data)
    set_sz <- colSums(gene.sets)
    
    unique_class_sz <- unique(set_sz)
    num_genes <- length(expr_data)
    
    divider <- c(1:num_genes)
    mean_table <- cumsum(expr_data)/divider
    mean_table_sq <- mean_table^2
    var_table <- cumsum(expr_data^2)/divider - (mean_table)^2 + pre.var
    class_sz_index <- rep(0,dim_sets[2])
    for (i in 1:dim_sets[2]){
        class_sz_index[i] <- which(set_sz[i]==unique_class_sz)
    }
    var_table <- var_table + pre.var
    hyge_stat <- .count_hyge_var_mean(dim_sets[1],unique_class_sz)
    
    z_var <- matrix(numeric(0),dim_sets[1],length(unique_class_sz))
    
    z_mean <- z_var
    max_value <- 1:dim_sets[1]
    for (j in 1:length(unique_class_sz)){
        prob_sum <- .count.prob.sum(dim_sets[1],hyge_stat$mean[,j],hyge_stat$var[,j])
        z_var[,j] <- 4*(var_table*prob_sum + mean_table_sq*hyge_stat$var[,j])
        z_mean[,j] <- mean_table*(2*hyge_stat$mean[,j]-max_value)}
    
    out = list(Z_var = z_var, Z_mean = z_mean, class_size_index = class_sz_index,var_table=var_table,mean_table_sq=mean_table_sq,set_sz=set_sz)
    return(out)
    
    
}

#########

.count.prob.sum <-
function(M, hyge.mean, hyge.var){
    tulos = matrix(rep(0,length(hyge.mean)),byrow=FALSE)
    N_tab = matrix(c(2:M),byrow=FALSE)
    tulos[1] = hyge.mean[1]
    tulos[2:M] = N_tab/(N_tab-1)*hyge.mean[2:M]-(hyge.mean[2:M]^2+hyge.var[2:M])/(N_tab-1)
    tulos = tulos
}

############

.count_hyge_var_mean <-
function(M,K){
    N = matrix(rep((1:M),length(K)),ncol=length(K))
    K = matrix(rep(K,M),byrow=TRUE,ncol=length(K))
    out1 = N*K/M
    out2 = N*(K*(M-K)/M^2)*((M-N)/(M-1))
    results = list(mean= out1, var= out2)
}

##############

.calc_z_var <-
function(num.genes,unique_class_sz_ln,pre_z_var,wgt2,var.constant){
    ones_matrix = matrix(1,num.genes,unique_class_sz_ln)
    ones_tmatrix = t(ones_matrix)
    median_matrix = apply(pre_z_var,2,median)
    pre_median_part = ones_tmatrix*median_matrix*wgt2
    median_part = t(pre_median_part)
    z_var = (pre_z_var + median_part + var.constant)^0.5
}

#############

.geneSetsList <-
function(data){
    data <- GSA.read.gmt(data)
    geneSets <- data$genesets
    names(geneSets) <- data$geneset.names
    return(geneSets)
}

###########

.des_perm4 <- function(expr.data, gene.sets, labels, rnaseq, comparisons=comparisons, perm.number){
    require(limma)
    require(Biobase)
    expr <- expr.data
    rownames(expr) <- NULL
    dime2 <- dim(expr)
    
    perms <- .perm4(labels,perm.number)
    perm.designs <- .getPointers(perms, labels)
    
    
    org.tscores <- list()
    ###### calculation of t-scores for original data ##########
    #design <- model.matrix(~0 + factor(labels))
    #colnames(design) <- paste("group",c(1:length(levels(factor(labels)))),sep="")
    design <- CreateDesignForLimma(colnames(expr.data),labels)
    groups <- colnames(design)
    
    
    
    if(rnaseq){
        print("RNAseq data inputed. This software expects raw counts!")
        dgelist <- DGEList(counts=expr,group=factor(labels)) # creating edgeR's data object
        dgelist <- calcNormFactors(dgelist, method="TMM") # TMM normalization
        y <- voom(counts=dgelist,design)
        fit1 <- lmFit(y, design)
    }
    
    else{
        
        mySet <- new("ExpressionSet",expr=as.matrix(expr))
        fit1 <- lmFit(mySet, design)
    }
    
    cont.matrix <- makeContrasts(contrasts=comparisons,levels=design)
    fit2 <- contrasts.fit(fit1, cont.matrix)
    fit2 <- eBayes(fit2)
    for(i in 1:length(comparisons)){
        org.tscores[[i]] <- fit2$t[,i]
    }
    
    ## Calculation of t-scores for permuted data for given comparisons ##
    perm.tscores <- vector("list",perm.number)
    for(k in 1:perm.number){
        perm.t <- list()
        design <- perm.designs[[k]]
        
        # generate permuted sample labels
        sample.labels.perm <- rep(NA,length(labels))
        grps <- c(1:length(unique(labels)))
        for(i in 1:length(grps)){
            sample.labels.perm[which(design[,i]==1)]= grps[i]
        }
        
        
        #colnames(design) <- paste("group",c(1:length(levels(factor(labels)))),sep="")
        colnames(design) <- groups
        
        if(rnaseq){
            dgelist <- DGEList(counts=expr,group=factor(sample.labels.perm)) # creating edgeR's data object
            dgelist <- calcNormFactors(dgelist, method="TMM") # TMM normalization
            y <- voom(counts=dgelist,design)
            fit1 <- lmFit(y, design)
        }
        
        
        else{
            mySet <- new("ExpressionSet",expr=as.matrix(expr))
            fit1 <- lmFit(mySet, design)
        }
        
        cont.matrix <- makeContrasts(contrasts=comparisons,levels=design)
        fit2 <- contrasts.fit(fit1, cont.matrix)
        fit2 <- eBayes(fit2)
        for(i in 1:length(comparisons)){
            perm.t[[i]] <- fit2$t[,i]
        }
        perm.tscores[[k]] <- perm.t
    }
    
    tscores.tot <- list()
    for(i in 1:length(comparisons)){
        tscores.tot[[i]] <- matrix(0,dime2[1],perm.number)
    }
    
    for(j in 1:length(comparisons)){
        for(k in 1:perm.number){
            tscores.tot[[j]][,k] <- perm.tscores[[k]][[j]]
        }
    }
    
    for(i in 1:length(comparisons)){
        tscores.tot[[i]] <- cbind(org.tscores[[i]], tscores.tot[[i]])
    }
    
    for(z in 1:length(comparisons)){
        rownames(tscores.tot[[z]]) <- rownames(expr.data)
        tscores.tot[[z]] <- aggregate(tscores.tot[[z]],by=list(rownames(tscores.tot[[z]])),FUN=mean)
        rownames(tscores.tot[[z]]) <- tscores.tot[[z]][,1]
        tscores.tot[[z]] <- tscores.tot[[z]][,2:dim(tscores.tot[[z]])[2]]
        tscores.tot[[z]][order(match(rownames(tscores.tot[[z]]),rownames(gene.sets))),]
    }
    
    return(tscores.tot)
    
    
}



############## perm4 ###########

.perm4 <- function(labels,perm.number){
    lbs <- list()
    u <- unique(labels)
    for(i in 1:(length(u))){
        lbs[[i]] <- which(labels==u[i])
    }
    perms <- vector("list",perm.number)
    for(j in 1:perm.number){
        lst <- list()
        lbs.org <- lbs
        lbs.ite <- lbs
        for(i in 1:length(lbs.org)){
            nonEmptyGrps <- sum(lapply(lbs.ite,function(x) length(x))!=0)
            if(nonEmptyGrps>=length(lbs.org[[i]])){
                lst[[i]] <- c(.resamp(unlist(lapply(lbs.ite,function(x) if(!length(x)==0){.resamp(x,1,replace=F)})),length(lbs.org[[i]]),replace=F))}
            else{
                lst[[i]] <- .resamp(unlist(lbs.ite),length(lbs.org[[i]]),replace=F)}
            lbs.ite <- .updateLbs(lbs.ite,lst[[i]])
        }
        perms[[j]] <- lst
    }
    
    return(perms)
}



## function to create design matrix for permuted data ##
.getPointers <- function(perms,labels){
    ngroups <- length(unique(labels))
    perm.labels <- list()
    d <- matrix(0,length(labels),ngroups)
    for(k in 1:length(perms)){
        for(i in 1:ngroups){
            n <- rep(0,length(labels))
            n[perms[[k]][[i]]] <- 1
            d[,i] <- n
        }
        perm.labels[[k]] <- d
    }
    return(perm.labels)
}

## function for sampling ##
.resamp <- function(x,...){if(length(x)==1) x else sample(x,...)}

##function to update sample labels ##
.updateLbs <- function(lbs,ind){
    l <- unlist(lbs)[-which(unlist(lbs)%in%ind)]
    lb <- list()
    for(i in 1:length(lbs)){
        lb[[i]] <- lbs[[i]][which(lbs[[i]]%in%l)]
    }
    lb
}


######################

.toMatrix <-
function(expr.data, gene.sets){
    if(is.matrix(gene.sets) | is.array(gene.sets) | is.data.frame(gene.sets) | is.table(gene.sets) ){
        gene.sets <- as.matrix(gene.sets)
    }
    
    if(is.list(gene.sets)){
        gene.set.matrix <- .listTOclMatrix(rownames(expr.data),gene.sets)
        gene.set.size <- dim(gene.set.matrix)
        colnames(gene.set.matrix) <- names(gene.sets)
        rownames(gene.set.matrix) <- rownames(expr.data)
        gene.sets <- gene.set.matrix
    }
    return(gene.sets)
}

#############

.listTOclMatrix <-
function(gene.vector, GO.list){
    num.classes <- length(GO.list)
    result <- list()
    for(i in 1:num.classes){
        result[[i]] <- gene.vector%in%GO.list[[i]]
    }
    out <- 1*(matrix(unlist(result),ncol=num.classes))
    return(out)
}

#############

.flipListStruct <- function(data){
    list_len <- length(data)
    max_val = 0
    ind_vals <- rep(0,list_len)
    loop_val <- 1
    for(k in 1:list_len){
        if(length(data[[k]])> 0){
            if(length(data[[k]])> 1 || data[[k]] > 0) {
                ind_vals[loop_val] <- k
                loop_val <- loop_val + 1
                max_val <- max(max_val, max(data[[k]]) )
            }
        }
        
    }
    ind_vals <- ind_vals[ind_vals > 0]
    out <- list(NULL)
    out[[max_val+1]] = 1
    out[[max_val+1]] = NULL
    for(k in ind_vals){
        for(l in data[[k]] ){
            out[[l]] <- c(out[[l]], k)
        }
    }
    out
}

#############
.emp <-
function(pos.data, neg.data){
    
    neg.dat.sz <- dim(neg.data)
    test_data <- rep(0, neg.dat.sz[2])
    for (k in 1:neg.dat.sz[2]){
        test_data[k] <- length(which(neg.data[,k] >= pos.data[k]))
    }
    test_data[test_data == 0] <- 0.5
    test_data <- (test_data)/(neg.dat.sz[1])
    #test_data <- -log10(test_data)
    test_data
    
    
}

###############

.KS.p.values <-
function(pos.data, perm.data){
    #library(MASS)
    tmp <- .mGSZ.adj.mean.std(pos.data, perm.data)
    pos.data <- tmp$pos.out
    perm.data <- tmp$perm.out
    rm(tmp)
    col.ind <- .pick.data.cols(perm.data)
    perm.data <- perm.data[,col.ind]
    pos.data <- pos.data[col.ind]
    emp.pval.class <- .emp(pos.data, perm.data)
    
    qvalue <- p.adjust(emp.pval.class,method="fdr")
    output <- list(EMP.class=emp.pval.class,class.ind=col.ind, FDR=qvalue)
    
    
    return(output)
}

###############

.logEVcdf <-
function(x, mu, sigma) {
    
    x <- (mu - x)/sigma  # Note! Flip the sign here
    out <- rep(0, length(x))
    if(!(is.vector(x)) ){
        d.sz <- dim(x)
        out <- matrix(out, d.sz[1], d.sz[2])
    }
    po1 <- which(x < 5) # switch for appr. vs. norm. eq.
    out[po1] <- -log(1 - exp(-exp(x[po1])) )
    x <- x[-po1]
    out[-po1] <- -x + exp(x)/2 - exp(2*x)/24 + exp(4*x)/2880 #Polyn. Expansion
    out <- out/log(10)
    out <- 10^(-out)
    out
}

###############

.logNORMcdf <-
function(x,mu,sigma){
    z <- (x-mu)/(sigma)
    
    #out <- rep(0,length(x))
    if(!(is.vector(x)) ){
        d.sz <- dim(x)
        out <- matrix(out, d.sz[1], d.sz[2])
    }
    out <- pnorm(-z,lower.tail=T)
    #out <- -log10(out)
    
    out
}

###############

.mAllez.p.values <-
function(pos.data, perm.data){
    #library(MASS)
    tmp <- .mGSZ.adj.mean.std(pos.data, perm.data)
    pos.data <- tmp$pos.out
    perm.data <- tmp$perm.out
    rm(tmp)
    col.ind <- .pick.data.cols(perm.data)
    perm.data <- perm.data[,col.ind]
    pos.data <- pos.data[col.ind]
    norm.param.all <- fitdistr(as.numeric(perm.data),'normal')$estimate
    norm.pval.class <- rep(0,length(pos.data))
    
    for (k in 1:length(pos.data)){
        norm.param.class <- fitdistr(as.vector(perm.data[,k]),'normal')$estimate
        norm.pval.class[k] <- .logNORMcdf(pos.data[k],norm.param.class[1],norm.param.class[2])
        
    }
    
    qvalue <- p.adjust(norm.pval.class,method="fdr")
    output <- list(NORM.class=norm.pval.class,class.ind=col.ind, FDR=qvalue)
    output
    
}

################

.mGSA.p.values <-
function(pos.data, perm.data) {
    #library(ismev)
    #library(MASS)
    tmp <- .mGSZ.adj.mean.std(pos.data, perm.data)
    scaled.pos.data <- tmp$pos.out
    scaled.perm.data <- tmp$perm.out
    rm(tmp)
    get.perm.p.vals = FALSE
    perm.dat.sz <- dim(scaled.perm.data)
    col.ind <- .pick.data.cols(scaled.perm.data) # discard cols with null var
    
    # First fit the parameters to the whole data
    scaled.perm.data <- scaled.perm.data[,col.ind]
    scaled.pos.data <- scaled.pos.data[col.ind]
    ev.p.val.class  <- rep(0, length(scaled.pos.data))
    
    for (k in 1:length(scaled.pos.data)){
        
        # Following fits parameters to perms of each class
        ev.param.class  <- gum.fit(as.vector(scaled.perm.data[,k]),show=FALSE)$mle
        ev.p.val.class[k]  <- .logEVcdf(scaled.pos.data[k], ev.param.class[1], ev.param.class[2]  )
        
    }
    qvalue <- p.adjust(ev.p.val.class,method="fdr")
    output <- list(EV.class = ev.p.val.class,class.ind=col.ind, FDR=qvalue)
    
    output
}

###############

.mGSZ.p.values <-
function(pos.data, perm.data) {
    #library(ismev)
    #library(MASS)
    tmp <- .mGSZ.adj.mean.std(pos.data, perm.data)
    pos.data <- tmp$pos.out
    perm.data <- tmp$perm.out
    rm(tmp)
    get.perm.p.vals = FALSE
    perm.dat.sz <- dim(perm.data)
    col.ind <- .pick.data.cols(perm.data)
    perm.data <- perm.data[,col.ind]
    pos.data <- pos.data[col.ind]
    ev.p.val.class  <- rep(0,length(pos.data))
    
    
    for (k in 1:length(pos.data)){
        
        # Following fits parameters to perms of each class
        ev.param.class  <- gum.fit(as.vector(perm.data[,k]),show=FALSE)$mle
        ev.p.val.class[k]  <- .logEVcdf(pos.data[k], ev.param.class[1], ev.param.class[2]  )
    }
    
    qvalue <- p.adjust(ev.p.val.class,method="fdr")
    output <- list(EV.class = ev.p.val.class,class.ind=col.ind,FDR=qvalue)
    
    
    output
}

###############

.SS.p.values <-
function(pos.data, perm.data){
    #library(MASS)
    tmp <- .mGSZ.adj.mean.std(pos.data, perm.data)
    
    scaled.pos.data <- tmp$pos.out
    scaled.perm.data <- tmp$perm.out
    rm(tmp)
    col.ind <- .pick.data.cols(scaled.perm.data)
    div.length <- length(scaled.perm.data)
    scaled.perm.data <- scaled.perm.data[,col.ind]
    scaled.pos.data <- scaled.pos.data[col.ind]
    emp.pval.class <- .emp(scaled.pos.data, scaled.perm.data)
    
    qvalue <- p.adjust(emp.pval.class,method="fdr")
    output <- list(EMP.class = emp.pval.class,class.ind=col.ind,FDR=qvalue)
    output
    
}

##############


.mGSZ.adj.mean.std <-
function(pos.data, perm.data) {
    perm.sz <- dim(perm.data)
    mean.prof <- apply(perm.data, 2, mean)
    std.prof    <- apply(perm.data, 2, sd)
    null.ind <- which(std.prof == 0)
    std.prof[null.ind] <- 0.1
    perm.out <- matrix(0, perm.sz[1], perm.sz[2])
    pos.out <- rep(0, perm.sz[2])
    for( k in 1:perm.sz[2]){
        perm.out[,k] <- ( perm.data[,k] - mean.prof[k]) /std.prof[k]
        pos.out[k]    <- ( pos.data[k] - mean.prof[k] )/std.prof[k]
    }
    out  <- list(pos.out = pos.out, perm.out = perm.out)
    out
}

##############

.pick.data.cols <-
function(data) {
    
    # Discard cols where var == 0
    # These disturb distribution fitting
    
    tmp <- apply(data, 2, var)
    out <- which(tmp > 0)
    out
}


####

.FC <- function(x, cl){
    x <- 2^x
    x.class1 <- x[(cl == 1)]
    x.class2 <- x[(cl == 2)]
    hoge <- log2(mean(x.class2)/mean(x.class1))
    return(hoge)
}

###############

.mGSZ.test.score <-
function(expr.data, gene.sets, wgt1, wgt2, pre.var, var.constant,start.val,...){
    
    num.genes <- length(expr.data)
    
    # Ordering of gene expression data and gene sets data
    
    gene.sets <- .toMatrix(expr.data,gene.sets)
    ord_out <- order(expr.data, decreasing= TRUE)
    expr.data <- expr.data[ord_out]
    set.dim <- dim(gene.sets)
    cols <- set.dim[2]
    gene.sets <- gene.sets[ord_out,]
    
    
    
    expr.data.ud <- expr.data[num.genes:1] # expression values turned up-side down
    # This does the analysis of the lower end
    
    
    num.genes=length(expr.data)
    num.classes=dim(gene.sets)[2]
    set_sz <- apply(gene.sets,2,sum)
    unique_class_sz_ln <- length(unique(set_sz))
    
    mAllez <- rep(0,num.classes)
    
    mgsa <- matrix(0,2,num.classes)
    
    pre_z_var.1 <- .sumVarMean_calc(expr.data, gene.sets, pre.var)
    pre_z_var.2 <- .sumVarMean_calc(expr.data.ud, gene.sets, pre.var)
    
    Z_var1 = .calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.1$Z_var,wgt2,var.constant)
    Z_var2 = .calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.2$Z_var,wgt2,var.constant)
    
    out2 = matrix(0,2,num.classes)
    mgsa1 <- numeric(num.classes)
    mgsa2 <- mgsa1
    mgsz1 <- mgsa1
    mgsz2 <- mgsa2
    
    
    for (k in 1:num.classes){
        po1 <- which(gene.sets[,k]==1)
        po0 <- which(gene.sets[,k]==0)
        if(length(po1) > 0 & length(po0) > 0){
            tmp1 = expr.data
            tmp1[po0] = 0
            tmp0 = expr.data
            tmp0[po1] = 0
            result1 = cumsum(tmp1)-cumsum(tmp0)-pre_z_var.1$Z_mean[,pre_z_var.1$class_size_index[k]]
            result2 = cumsum(tmp1[num.genes:1])-cumsum(tmp0[num.genes:1])-pre_z_var.2$Z_mean[,pre_z_var.2$class_size_index[k]]
            
            result1[1:start.val] <- 0
            result2[1:start.val] <- 0
            A = result1/Z_var1[,pre_z_var.1$class_size_index[k]]
            B = result2/Z_var2[,pre_z_var.2$class_size_index[k]]
            
            mAllez[k] <- abs(A[num.genes])
            mgsa1[k] <- A[round(num.genes/2)]
            mgsa2[k] <- B[round(num.genes/2)]
            mgsz1[k] <- max(abs(A))
            mgsz2[k] <- max(abs(B))
        }
    }
    mgsa[1,]<-mgsa1
    mgsa[2,]<-mgsa2
    out2[1,]<-mgsz1
    out2[2,]<-mgsz2
    
    
    d.ind <- apply(out2,2,which.max) # index of the max value (1:up and 2:down)
    gsz.direction <- d.ind
    gsz.direction[gsz.direction==1]="up"
    gsz.direction[gsz.direction==2]="down"
    
    result1 = list(mGSZ.scores = apply(out2,2,max), mGSZ.direction=gsz.direction, upper.tail = mgsz1, lower.tail=mgsz2, mAllez.scores = mAllez, mGSA.scores= apply(abs(mgsa),2,max))
    result2 = list(Z_var1=Z_var1,Z_var2=Z_var2,Z_mean1=pre_z_var.1$Z_mean,Z_mean2=pre_z_var.2$Z_mean,class_size_index1=pre_z_var.1$class_size_index,class_size_index2=pre_z_var.2$class_size_index)
    
    out = list(gene.set.scores=result1, var.attributes=result2)
    return(out)
}

#######
.trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#######

.KS.score <-
function(expr.data, gene.sets){
    ord <- order(expr.data,decreasing=FALSE)
    expr.data <- expr.data[ord]
    gene.sets <- gene.sets[ord,]
    set.dim <- dim(gene.sets)
    result1 <- rep(0,set.dim[2])
    result2 <- result1
    profile2 <- matrix(0,set.dim[2],set.dim[1])
    expr.data.len <- length(expr.data)
    for(i in 1:set.dim[2]){
        po1 <- which(gene.sets[,i]==1)
        po0 <- which(gene.sets[,i]==0)
        tmp1 <- expr.data
        tmp1[po0] <- 0
        tmp0 <- expr.data
        tmp0[po1] <- 0
        tmp2 <- cumsum(abs(tmp1))/sum(abs(tmp1))
        tmp2b <- tmp1
        tmp2b[po1] <- 1
        tmp2c <- cumsum(tmp2b)/sum(tmp2b)
        tmp0b <- tmp0
        tmp0b[po0] <- 1
        tmp3  <- cumsum(abs(tmp0))/sum(abs(tmp0))
        tmp3b <- cumsum(abs(tmp0b))/sum(abs(tmp0b))
        result1[i]  = max(abs(tmp2c - tmp3b))
        result2[i] = max(abs(tmp2 - tmp3b))
    }
    
    out <- list(KS.org = result1,KS.mod = result2)
    
}
################

.sumTestscore <-
function(expr.data, gene.sets){
    num.classes <- dim(gene.sets)[2]
    Mu <- mean(expr.data)
    sum.score = rep(0,num.classes)
    for(i in 1:num.classes){
        pos1 <- which(gene.sets[,i]==1)
        ##Calculation of sum of squares
        sum.score[i] <- abs(sum(expr.data[pos1]-Mu))
        #sum.score[i] <- sum(expr.data[pos1])
        
    }
    return(sum.score)
}

##############

.SS.test.score <-
function(expr.data, gene.sets){
    num.classes <- dim(gene.sets)[2]
    Mu <- mean(expr.data)
    SS.score = rep(0,num.classes)
    for(i in 1:num.classes){
        pos1 <- which(gene.sets[,i]==1)
        ##Calculation of sum of squares
        SS.score[i] <- sum((expr.data[pos1]-Mu)^2)
        #SS.score[i] <- sum((expr.data[pos1])^2)
        
    }
    return(SS.score)
}

##########

CreateDesignForLimma <- function( ColNames, DesignNames){
    #ColNames = names of columns in expression data
    #DesignNames = corresponding sample names
    UniqNames = unique( DesignNames)
    Output = matrix(0, length( ColNames), length(UniqNames))
    for( k in 1:length(UniqNames)){
        Ind = which(DesignNames == UniqNames[k])
        Output[Ind,k] = 1
    }
    colnames( Output) = UniqNames
    rownames( Output) = ColNames
    Output
}

#############


