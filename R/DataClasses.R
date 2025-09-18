

setClassUnion("matrixORdata.frame", c("matrix", "data.frame"))
setClassUnion("matrixORdata.frameORlist", c("matrix", "data.frame","list"))
setClassUnion("numericORcharacter", c("numeric", "character"))

GSZ <- setClass("GSZ",slots=c(exprdat="matrixORdata.frame",genesets="matrixORdata.frameORlist",sampleLabels="numericORcharacter"))


mGSZResults <-
    setClass("mGSZResult",
            contains="GSZ",
            slots=c(mGSZ ="data.frame",sampleLabels="numericORcharacter",exprdat="matrixORdata.frame",genesets="matrixORdata.frameORlist"))

mGSZResultsAll <-
setClass("mGSZResultAll",
contains="GSZ",
slots=c(mGSZ ="data.frame",mGSA = "data.frame", mAllez = "data.frame",wKS = "data.frame",SUM = "data.frame",SS ="data.frame",sampleLabels="numericORcharacter",exprdat="matrixORdata.frame",genesets="matrixORdata.frameORlist"))



mGSZmResults <-
setClass("mGSZmResult",
contains="GSZ",
slots=c(results="list",AllComparisons="character"))

