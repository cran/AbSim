#' Clonally expands the simulated repertoire generated from fullRepertoire function.
#' @param ab.repertoire The output from the fullRepertoire function. This takes should
#' be a nested list, with the first element containing a list of sequence arrays,
#' the second element containing the corresponding list of names, and the third
#' element containing a list of trees. The following list elements (containing the
#' sampling information) will be disregarded for repertoire expansion.
#' @param rep.size Controls the total size of the final repertoire. This number
#' should be greater than the total number of sequences that are provided as input.
#' Note that currently using the "powerlaw" distribution will return less than the exact
#' rep.size parameter due to integer values of clonal frequencies.
#' @param distribution This parameter controls how the clonal frequency of the
#' repertoire is distributed. Options include "powerlaw" ("pl"), "identical" ("id")
#' @return Returns a list with three elements. The first list element contains
#' a character array, with the sequences composing the repertoire. The second element
#' is a character array with the names of the sequences, and the third element in the list
#' corresponds to the original phylogenetic trees that served as the basis for expansion
#' @param with.germline Logical - If false, the germline from each lineage will be removed.
#' @param dist.parameters Supplies the parameters for how the repertoire should
#' be distributed.
#' @seealso fullRepertoire
#' @export
#' @examples
#' \dontrun{
#' clonalExpansion(ab.repertoire=fullRepertoire.output,
#'  rep.size=3*length(unlist(fullRepertoire.output[[1]])),
#'  distribution="identical",
#'  with.germline="FALSE")
#'}

clonalExpansion <- function(ab.repertoire,
                            rep.size,
                            distribution,
                            with.germline,
                            dist.parameters){
  if(missing(dist.parameters)==TRUE) dist.parameters <- NULL
  output_list <- list()
  if (with.germline==FALSE){
    for(i in 1:length(ab.repertoire[[1]])){
      ab.repertoire[[1]][[i]] <- ab.repertoire[[1]][[i]][-1]
      ab.repertoire[[2]][[i]] <- ab.repertoire[[2]][[i]][-1]
    }
  }
  output_list[[1]] <- unlist(ab.repertoire[[1]])
  output_list[[2]] <- unlist(ab.repertoire[[2]])
  output_list[[3]] <- ab.repertoire[[3]]
  input_length <- length(output_list[[1]])
  if(length(unlist(ab.repertoire[[1]]))>=rep.size){
    print("Warning: number of core sequences is greater than specified repertoire size")
    return(output_list)
  }
  else{
    if(distribution=="id" || distribution=="identical"){
      clone_freq <- as.integer(rep.size/input_length)
      output_list[[1]] <- rep(output_list[[1]],clone_freq)
      output_list[[2]] <- rep(output_list[[2]],clone_freq)
    }
    else if(distribution=="powerlaw" || distribution=="pl"){
      if(is.null(dist.parameters)==TRUE){
        stop("Supply alpha value in dist.parameters argument")

      }
      xmin <- 1
      x <- input_length
      alpha <- dist.parameters
      dist <- poweRlaw::dpldis(1:x, xmin, alpha)
      print(dist)
      new_seq_array <- c()
      new_name_array <- c()
      print(length(output_list[[1]]))
      for(i in 1:length(output_list[[1]])){
        new_seq_array <- append(new_seq_array,rep(output_list[[1]][i],as.integer(dist[i]*rep.size)))
        new_name_array <- append(new_name_array,rep(output_list[[2]][i],as.integer(dist[i]*rep.size)))
      }
      output_list[[1]] <- new_seq_array
      output_list[[2]] <- new_name_array
    }

  }
  return(output_list)

}
