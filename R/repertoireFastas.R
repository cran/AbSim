#' Converts AbSim's output sequences and names into fasta file format. Each tree
#' within the repertire will be stored in an individual fasta file. Utilizes
#' sink() and cloneAllConnections() from R base
#' @param ab.repertoire The output from the either the singleLineage function or
#' the fullRepertoire function. This takes should be a nested list, with the first
#' element containing a list of sequence arrays, the second element containing the
#' corresponding list of names, and the third element containing a list of trees.
#' The following list elements (containing the sampling information) will be
#' disregarded for repertoire expansion.
#' @param fasta.name The name of the fasta files. Each tree will be numbered
#' corresponding to the index within the ab.repertoire list. For example,
#' the second tree in an ab.repertoire object will be "fasta.name_2.fasta",
#' @param with.germline Logical value specifying if the germline sequences
#' should be included in the fasta file (as the first sequence)
#' @export
#' @examples
#' \dontrun{
#' repertoireFastas(ab.repertoire=ab.repertoire.output,
#'  fasta.name="sample_trees",
#'  with.germline=TRUE
#'}

repertoireFastas <- function(ab.repertoire,
                             directory.string,
                             fasta.name,
                             with.germline){
  for(i in 1:length(ab.repertoire[[1]])){
    ape::write.dna(x=as.matrix(x=c(ab.repertoire[[1]][[i]],
                                  ab.repertoire[[2]][[i]])),
                   format="fasta",
                   file=paste(directory.string,"/",fasta.name,"_",i,".fasta",sep=""))
  }
}

repertoireFastas <- function(ab.repertoire,
                             fasta.name,
                             with.germline){
  init.count <- 1
  for(i in 1:length(ab.repertoire[[1]])){
    unlisted.rep.seq <- unlist(ab.repertoire[[1]][[i]])
    unlisted.rep.names <- unlist(ab.repertoire[[2]][[i]])
    sink(paste(fasta.name,"_",i,".fasta",sep=""))
    if(with.germline==FALSE) init.count <- 2
    for(j in init.count:length(unlisted.rep.seq)){
      cat(paste(">",unlisted.rep.names[j],sep=""))
      cat("\n")
      cat(unlisted.rep.seq[j])
      if(j!=length(unlisted.rep.seq)) cat("\n")

    }
    closeAllConnections()
  }
}


