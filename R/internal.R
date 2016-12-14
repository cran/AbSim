#internal functions

.VDJ_RECOMBIN_FUNCTION <- function(v_seq, d_seq, j_seq, method){
  base_array <- c("a", "t", "g", "c")
  base_probs <- c(.25,.25,.25,.25)
  insertion_function <- function(seq1, seq2){

  }
  vdj_deletion_function <- function(seq1,seq2){

  }
  ## need to recombine V and D and D and J
  if (method=="naive"){
    base_array <- c("a", "t", "g", "c")
    base_probs <- c(.25,.25,.25,.25)
    vdj_options <- c("none", "deletion", "insertion")
    vdj_options_prob <- c(.4,.3,.3)
    recomb_decision <- sample(x = vdj_options, 2, replace=TRUE, prob = vdj_options_prob)
    while(recomb_decision[1]=="none" && recomb_decision[2]=="none"){
      recomb_decision <- sample(x = vdj_options, 2, replace=TRUE, prob = vdj_options_prob)
    }
    ############ starts the combination between V and D
    if(recomb_decision[1]=="none"){ # no change
      v_d <- paste(v_seq, d_seq, sep="")
    }
    else if(recomb_decision[1]=="deletion"){
      # v_seq_new loses 3' end and d_seq_new1 loses up to 5' bases
      v_seq_new <- substring(v_seq,first=1, last=nchar(v_seq)-sample(x=c(0,1,2,3,4,5),1,replace=TRUE, prob=c(0.3,0.2,0.2,0.1,0.1,0.1)))
      d_seq_new5 <- substring(d_seq, first=sample(x=c(1,2,3,4,5,6),1,replace=TRUE,prob=c(0.3,0.2,0.2,0.1,0.1,0.1)), last=nchar(d_seq))
      v_d <- paste(v_seq_new, d_seq_new5, sep="")
    }
    else if(recomb_decision[1]=="insertion"){
      insertion_length <- sample(x=c(1,2,3,4,5),1,replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2))
      insertion_array <- sample(x=base_array, insertion_length, replace=TRUE, base_probs)
      insertion_string <- paste(insertion_array, collapse='')
      v_d <- paste(v_seq, insertion_string, d_seq, sep="")

    } ### END recombination between v and d regions and start D and J
    if(recomb_decision[2]=="none"){
      vdj <- paste(v_d, j_seq, sep="")
    }
    else if(recomb_decision[2]=="deletion"){
      ### need to cut off the end of v_d and cut off the 5' of J
      v_d_new <- substring(v_d,first=1, last=nchar(v_d)-sample(x=c(1,2,3,4,5),1,replace=TRUE, prob=c(0.2,0.2,0.2,0.2,0.2)))
      j_new <- substring(j_seq, first=sample(x=c(1,2,3,4,5),1,replace=TRUE,prob=c(0.2,0.2,0.2,0.2,0.2)), last=nchar(j_seq))
      vdj <- paste(v_d_new, j_new, sep="")
    }
    else if(recomb_decision[2]=="insertion"){
      dj_insertion_length <- sample(x=c(3,6,9,12),1,replace=TRUE, prob=c(0.25,0.25,0.25,0.25))
      dj_insertion_array <- sample(x=base_array, dj_insertion_length, replace=TRUE, base_probs)
      dj_insertion_string <- paste(dj_insertion_array, collapse='')
      vdj <- paste(v_d, dj_insertion_string, j_seq, sep="")
    } ### ends the second recom_decision if-else


    return(vdj)
  }
  ## use data driven method
}

.SHM_FUNCTION_SEQUENCE4 <- function(vdj_seq, mut_param,v_seq,
                                   d_seq,j_seq, SHM.nuc.prob){
  #hot_spot_df <- read.table(file="~/PHD/Phylo/simu_results/pipeline/fivemer_one.csv", header=TRUE)
  #hot_spot_df <- read.table(file="~/PHD/Phylo/simu_results/pipeline/fivemer.csv", header=TRUE)
  #one_spot_df <- read.table(file="~/PHD/Phylo/simu_results/pipeline/fivemer_one.csv", header=TRUE)
  # @vdj_seq = sequence that will be mutated
  # @mut_param = determines how the mutations are given to the sequence
  # "naive" for random in whole VDJ, "cdr" for more in cdr region, "hotspot"
  ## 10^-3/bp/generation -

  base_line_mutations <- 0
  if(mut_param=="naive" || mut_param=="all"){
    holding_mut <- sample(x=c(0,1), nchar(vdj_seq), replace=TRUE,
                          c(SHM.nuc.prob[1], 1-SHM.nuc.prob[1]))
    for (i in 1:nchar(vdj_seq)){
      if(holding_mut[i]==0){
        holding_char <- substr(vdj_seq, i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(rep(1/3,3)))
        #substr(vdj_seq, i, i) <- sample(x=c("A", "T", "G","C"), 1, replace=TRUE, c(rep(.25,4)))
        base_line_mutations <- base_line_mutations + 1
      }
    }
  }
  if(mut_param=="data" || mut_param=="all"){
    index_CDR_start <- nchar(v_seq)-15
    index_CDR_stop <- nchar(vdj_seq)
    CDR_length <- index_CDR_stop - index_CDR_start
    if(mut_param=="data") CDR_prob <- SHM.nuc.prob
    else if(mut_param=="all") CDR_prob <- SHM.nuc.prob[2]
    no_CDR_prob <- 1-CDR_prob
    #CDR_mut <- sample(x=c(0,1), nchar(CDR_length), replace=TRUE, c(no_CDR_prob,CDR_prob))
    for(i in 81:114){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
      }
    }
    for(i in 168:195){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
      }
    }
    for(i in index_CDR_start:index_CDR_stop){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
      }
    }
  }
  if(mut_param == "motif" || mut_param == "all"){
    random_spot <- hotspot_df[sample(nrow(hotspot_df)),]
    current_hot_spots <- 0
    if(mut_param=="motif") hot_spot_limit <- SHM.nuc.prob
    else if(mut_param=="all") hot_spot_limit <- SHM.nuc.prob[3]
    #hot_spot_limit <- SHM.nuc.prob
    for(i in 1:nrow(random_spot)){
      if(grepl(pattern = random_spot$pattern[i],x = vdj_seq)){
        current_hot_spots <- current_hot_spots+1
      }
      if(current_hot_spots>hot_spot_limit) break
      vdj_seq <- gsub(random_spot$pattern[i],replacement = paste(substring(random_spot$pattern[i],first=1,last=2),sample(x = c("A","C","G","T"),1,replace=TRUE,prob = c(100,0,0,0)),substring(random_spot$pattern[i],first=4,last=5),sep=""),x = vdj_seq)

    }
  }
  if(mut_param == "wrc" || mut_param == "all"){
    random_spot <- one_spot_df[sample(nrow(one_spot_df)),]
    current_hot_spots <- 0
    if(mut_param=="motif") hot_spot_limit <- SHM.nuc.prob
    else if(mut_param=="all") hot_spot_limit <- SHM.nuc.prob[4]
    #hot_spot_limit <- SHM.nuc.prob
    for(i in 1:nrow(random_spot)){
      if(grepl(pattern = random_spot$pattern[i],x = vdj_seq)){
        current_hot_spots <- current_hot_spots+1
      }
      if(current_hot_spots>hot_spot_limit) break
      vdj_seq <- gsub(random_spot$pattern[i],replacement = paste(substring(random_spot$pattern[i],first=1,last=2),sample(x = c("A","C","G","T"),1,replace=TRUE,prob = c(100,0,0,0)),substring(random_spot$pattern[i],first=4,last=5),sep=""),x = vdj_seq)
    }
  }
  return(vdj_seq)

}


.applyBaseLine <- function(sequence_array, mutation_prob){
  no_mutation <- 1-mutation_prob
  for(i in 1:length(sequence_array)){
    holding_base <- sample(x=c(0,1), nchar(sequence_array[i]),replace=TRUE,
                           c(no_mutation, mutation_prob))
    for(j in 1:nchar(sequence_array[i])){
      if(holding_base[j]==1){
        holding_char <- substr(sequence_array[i],j,j)
        if(holding_char=="A") substr(sequence_array[i], j, j) <- sample(x=c("T", "G","C"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="T") substr(sequence_array[i], j, j) <- sample(x=c("A", "G","C"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="G") substr(sequence_array[i], j, j) <- sample(x=c("A", "T","C"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="C") substr(sequence_array[i], j, j) <- sample(x=c("A", "G","T"), 1, replace=TRUE, c(rep(1/3,3)))
        #substr(sequence_array[i], j, j) <- sample(x=c("A", "T", "G","C"), 1, replace=TRUE, c(rep(.25,4)))
      }
    }
  }
  return(sequence_array)
}


.branchingProcess3 <- function(input_string, branch_node, new_node, node_type){
  if(node_type=="lineage"){
    if (grepl(pattern = paste(",",as.character(branch_node),")",sep=""),x=input_string)==TRUE){
      output_string <- gsub(paste(",",as.character(branch_node),")",sep=""), paste(",","(", as.character(branch_node),",",as.character(new_node),"))",sep = ""),input_string)
    }

    else if(grepl(pattern = paste("(",(branch_node),",",sep=""),x=input_string,fixed=TRUE)==TRUE){
      output_string <- gsub(paste("(",as.character(branch_node),",",sep=""), paste("((", as.character(branch_node),",",as.character(new_node),"),",sep = ""),input_string,fixed=TRUE)

    }
  }
  else if(node_type=="SHM"){
    if (grepl(pattern = paste(",",as.character(branch_node),")",sep=""),x=input_string)==TRUE){
      output_string <- gsub(paste(",",as.character(branch_node),")",sep=""), paste(",","(", as.character(branch_node),",",as.character(new_node),"))",sep = ""),input_string)
    }

    else if(grepl(pattern = paste("(",(branch_node),",",sep=""),x=input_string,fixed=TRUE)==TRUE){
      output_string <- gsub(paste("(",as.character(branch_node),",",sep=""), paste("((", as.character(branch_node),",",as.character(new_node),"),",sep = ""),input_string,fixed=TRUE)

    }
  }
  return(output_string)
}
