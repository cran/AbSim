#' Simulates full heavy chain antibody repertoires for either human or mice.
#' @param max.tree.num Integer value describing maximum number of trees allowed
#'  to generate the core sequences of the repertoire. Each of these trees is started
#'  by an independent VDJ recombination event.
#' @inheritParams singleLineage
#' @return Returns a nested list. output[[1]][[1]] is an array of the simulated sequences
#' output[[2]][[1]] is an array names corresponding to each sequence. For example, output[[2]][[1]][1]
#' is the name of the sequence corresponding to output[[1]][[1]][1]. The simulated tree of this is found in
#' output[[3]][[1]]. The length of the output list is determined by the number of sampling points
#' Thus if you have two sampling points, output[[4]][[1]] would be a character array holding the sequences
#' with output[[5]][[1]] as a character array holding the corresponding names. Then the sequences recovered
#' second sampling point would be stored at output[[6]][[1]], with the names at output[[7]][[1]]. This
#' nested list was designed for full antibody repertoire simulations, and thus, may seem unintuitive
#' for the single lineage function. The first sequence and name corresponds to the germline sequence
#' that served as the root of the tree. See vignette for comprehensive example
#' @export
#' @seealso singleLineage
#' @examples
#' fullRepertoire(max.seq.num=51,max.timer=150,
#'  SHM.method="naive",baseline.mut = 0.0008,
#'  SHM.branch.prob = "identical", SHM.branch.param = 0.05,
#'  SHM.nuc.prob = 15/350,species="mus",
#'  VDJ.branch.prob = 0.1,proportion.sampled = 1,
#'  sample.time = 50,max.tree.num=3)

fullRepertoire <- function(max.seq.num,
                          max.timer,
                          SHM.method,
                          baseline.mut,
                          SHM.branch.prob,
                          SHM.branch.param,
                          SHM.nuc.prob,
                          species,
                          VDJ.branch.prob,
                          proportion.sampled,
                          sample.time,
                          max.tree.num
                          ){
  tree_list <- c()
  tree_list[1:max.tree.num] <- "(0,1)L;"
  seq_list <- list()
  name_list <- list()

  seq_per_tree <- rep(1,max.tree.num)
  germline_vseq <- c()
  germline_dseq <- c()
  germline_jseq <- c()
  germline_seq <- c()
  germline_name <- c()
  #tree_list <- list()
  max.seq.num <- max.seq.num + 1
  output_list <- list()
  max_SHM <- max.seq.num - max.tree.num
  if(species=="mus" || species=="mouse"){
    indV <- sample(x = 1:nrow(blc6_v_df),size = 1,replace = FALSE)
    indD <- sample(x = 1:nrow(blc6_d_df),size=1,replace=FALSE)
    indJ <- sample(x = 1:nrow(blc6_j_df),size=1,replace=FALSE)
    germline_name[1] <- paste(as.character(blc6_v_df[[1]][indV]),as.character(blc6_d_df[[1]][indD]),as.character(blc6_j_df[[1]][indJ]),sep="")
    germline_vseq[1] <- as.character(blc6_v_df[[2]][indV])
    germline_dseq[1] <- as.character(blc6_d_df[[2]][indD])
    germline_jseq[1] <- as.character(blc6_j_df[[2]][indJ])
    germline_seq[1] <- paste(as.character(blc6_v_df[[2]][indV]),as.character(blc6_d_df[[2]][indD]),as.character(blc6_j_df[[2]][indJ]),sep="")
    seq_list[[1]] <- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(blc6_v_df[[2]][indV]), as.character(blc6_d_df[[2]][indD]),as.character(blc6_j_df[[2]][indJ]),"naive"))
    name_list[[1]] <- paste(paste("S","1", sep=""),1,sep="_")
  }
  else if(species=="hum" || species== "human"){
    indV <- sample(x = 1:nrow(hum_v_df),size = 1,replace = FALSE)
    indD <- sample(x = 1:nrow(hum_d_df),size=1,replace=FALSE)
    indJ <- sample(x = 1:nrow(hum_j_df),size=1,replace=FALSE)
    germline_name[1] <- paste(as.character(hum_v_df[[1]][indV]),as.character(hum_d_df[[1]][indD]),as.character(hum_j_df[[1]][indJ]),sep="")
    germline_vseq[1] <- as.character(hum_v_df[[2]][indV])
    germline_dseq[1] <- as.character(hum_d_df[[2]][indD])
    germline_jseq[1] <- as.character(hum_j_df[[2]][indJ])
    germline_seq[1] <- paste(as.character(hum_v_df[[2]][indV]),as.character(hum_d_df[[2]][indD]),as.character(hum_j_df[[2]][indJ]),sep="")
    seq_list[[1]] <- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(hum_v_df[[2]][indV]), as.character(hum_d_df[[2]][indD]),as.character(hum_j_df[[2]][indJ]),"naive"))
    name_list[[1]] <- paste(paste("S","1", sep=""),1,sep="_")
  }

  current_seq_count <- 1
  VDJ_count <- 1
  SHM_count <- 0
  next_node <- 2
  output_list <- list()
  if(class(sample.time)=="numeric" && length(sample.time)>0){
    for(i in seq(1,length(sample.time),2)){
      output_list[[i+3]] <- list()
      output_list[[i+4]] <- list()
    }
  }

  seq_type <- c()
  new_VDJ_prob <- VDJ.branch.prob
  if(SHM.branch.prob=="identical"){
    new_SHM_prob <- rep(SHM.branch.param,max.seq.num)
  }
  else if(SHM.branch.prob=="uniform"){
    new_SHM_prob <- stats::runif(n=max.seq.num,min=SHM.branch.param[1],max=SHM.branch.param[2])
  }
  else if(SHM.branch.prob=="exponential"||"exp"){
    new_SHM_prob <- stats::rexp(n=max.seq.num,rate=SHM.branch.param)
  }
  else if(SHM.branch.prob=="lognorm"|| "lognormal"){
    new_SHM_prob <- stats::rlnorm(n=max.seq.num,meanlog = SHM.branch.param[1],sdlog = SHM.branch.param[2])
  }
  else if(SHM.branch.prob=="normal"||"norm"){
    new_SHM_prob <- stats::rnorm(n=max.seq.num,mean=SHM.branch.param[1],sd=SHM.branch.param[2])
  }

  sample_seq_list <- list()
  sample_names_list <- list()
  output_tree_text <- "(0,1)L;"
  sample_index <- 1
  sample.time <- sort(sample.time,decreasing = FALSE)
  current_tree_num <- 1
  for(i in 1:max.timer){
    if(length(unlist(seq_list))>=max.seq.num) break
    if(i==sample.time[sample_index] && length(unlist(seq_list)) >= 1){
      for(z in 1:length(seq_list)){
        current_len <- length(seq_list[[z]])
        holding_size <- as.integer(proportion.sampled * current_len)
        holding_ind <- sample.int(n = current_len,
                                  size=holding_size,
                                  replace = FALSE,
                                  prob = rep(x = 1/current_len,current_len))
        output_list[[2*sample_index+2]][[z]] <- seq_list[[z]][holding_ind]
        output_list[[2*sample_index+3]][[z]] <- paste(name_list[[z]][holding_ind],sample.time[sample_index],sep="_")
      }
    if(length(sample.time) > sample_index) sample_index <- sample_index + 1
    }
    if(length(unlist(seq_list))>=2){
      for(u in 1:length(seq_list)){
        seq_list[[u]] <- .applyBaseLine(seq_list[[u]],baseline.mut)
      }
    }
    is_new_VDJ <- sample(x=c(0,1), replace=TRUE,size = 1, prob=c(VDJ.branch.prob,
                                                        1- VDJ.branch.prob))
    if(is_new_VDJ==0 && current_seq_count<max.seq.num && current_tree_num<max.tree.num){
      if(species=="mus" || species=="mouse"){
        indV <- sample(x = 1:nrow(blc6_v_df),size = 1,replace = FALSE)
        indD <- sample(x = 1:nrow(blc6_d_df),size=1,replace=FALSE)
        indJ <- sample(x = 1:nrow(blc6_j_df),size=1,replace=FALSE)
        germline_name[current_tree_num+1] <- paste(as.character(blc6_v_df[[1]][indV]),as.character(blc6_d_df[[1]][indD]),as.character(blc6_j_df[[1]][indJ]),sep="")
        germline_seq[current_tree_num+1] <- paste(as.character(blc6_v_df[[2]][indV]),as.character(blc6_d_df[[2]][indD]),as.character(blc6_j_df[[2]][indJ]),sep="")
        seq_list[[current_tree_num+1]] <- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(blc6_v_df[[2]][indV]), as.character(blc6_d_df[[2]][indD]),as.character(blc6_j_df[[2]][indJ]),"naive"))
      }
      else if(species=="hum" || species== "human"){
        indV <- sample(x = 1:nrow(hum_v_df),size = 1,replace = FALSE)
        indD <- sample(x = 1:nrow(hum_d_df),size=1,replace=FALSE)
        indJ <- sample(x = 1:nrow(hum_j_df),size=1,replace=FALSE)
        germline_name[current_tree_num+1] <- paste(as.character(hum_v_df[[1]][indV]),as.character(hum_d_df[[1]][indD]),as.character(hum_j_df[[1]][indJ]),sep="")
        germline_seq[current_tree_num+1] <- paste(as.character(hum_v_df[[2]][indV]),as.character(hum_d_df[[2]][indD]),as.character(hum_j_df[[2]][indJ]),sep="")
        seq_list[[current_tree_num+1]] <- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(hum_v_df[[2]][indV]), as.character(hum_d_df[[2]][indD]),as.character(hum_j_df[[2]][indJ]),"naive"))
      }
      name_list[[current_tree_num+1]] <- paste(paste("S",next_node, sep=""),i,sep="_")

      current_seq_count <- current_seq_count + 1
      current_tree_num <- current_tree_num + 1
      next_node <- next_node + 1
      seq_per_tree[current_tree_num+1] <- seq_per_tree[current_tree_num+1]
      tree_list[current_tree_num] <- gsub(pattern="1",replacement = as.character(next_node-1),x=tree_list[current_tree_num])
    }
    if(next_node>= max.seq.num) break
    if(current_seq_count>=1){
      for(o in 1:length(seq_list)){
       for(j in 1:length(seq_list[[o]])){
        is_new_SHM <- sample(x=c(0,1),size=1, replace=TRUE, prob=c(new_SHM_prob[next_node], 1- new_SHM_prob[next_node]))
        if (is_new_SHM==0 && next_node<max.seq.num && SHM_count<max_SHM ){
          holding_jsub <- gsub("_.*","",name_list[[o]][j])
          holding_jsub <- gsub("[^0-9]","",holding_jsub)
          seq_list[[o]][seq_per_tree[o]+1] <- .SHM_FUNCTION_SEQUENCE4(seq_list[[o]][j],
                                                          SHM.method,germline_vseq[o],germline_dseq[o],germline_jseq[o],SHM.nuc.prob)
          name_list[[o]][seq_per_tree[o]+1] <- paste(paste("S",next_node, sep=""),i,sep="_")
          tree_list[o] <- .branchingProcess3(tree_list[o],holding_jsub,next_node,"SHM")
          next_node <- next_node + 1
          current_seq_count <- current_seq_count + 1
          SHM_count <- SHM_count + 1
          seq_per_tree[o] <- seq_per_tree[o] + 1

        }
        holding_no_nas <- unlist(seq_list)
        if(length(holding_no_nas[is.na(holding_no_nas)==FALSE])>=max.seq.num) break
       }
      holding_no_nas <- unlist(seq_list)
      if(length(holding_no_nas[is.na(holding_no_nas)==FALSE])>=max.seq.num) break      }
    }
    holding_no_nas <- unlist(seq_list)
    if(length(holding_no_nas[is.na(holding_no_nas)==FALSE])>=max.seq.num) break
  }
  for(i1 in 1:length(seq_list)){
    output_list[[1]][[i1]] <- append(germline_seq[i1],seq_list[[i1]])
    output_list[[2]][[i1]] <- append(germline_name[i1],name_list[[i1]])
    temp.tree <- ape::read.tree(text=tree_list[i1])
    output_list[[3]][[i1]] <- temp.tree
    for(i2 in 2:length(output_list[[3]][[i1]]$tip.label)){
      output_list[[3]][[i1]]$tip.label[1] <- output_list[[2]][[i1]][1]
      output_list[[3]][[i1]]$tip.label[i2] <- output_list[[2]][[i1]][grep(output_list[[2]][[i1]],pattern=paste("S",output_list[[3]][[i1]]$tip.label[i2],"_",sep=""))]
    }
   }

  return(output_list)
}
