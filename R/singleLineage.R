#' Antibody lineage simulations using only one set of V(D)J germline genes. The main difference between this function
#' and the fullRepertoire function is that there can be multiple VDJ recombination events within one tree. Each VDJ recombination
#' event in the singleLineage function is a branching event within the existing tree, whereas the VDJ recombination
#' events in the fullRepertoire function start a new tree.
#'
#' @param max.seq.num The maximum number of tips allowed at the end of the simulation. The simulation
#' will end when either this or the max.timer is reached. Note - this function does not take clonal
#' frequency into account. This parameter resembles the species richness, or the measure of unique sequences
#' in the repertoire.
#' @param max.timer The maximum number of time steps allowed during the simulation. The simulation
#' will end when either this or the max.seq.num is reached.
#' @param SHM.method The mode of SHM speciation events. Options are either: "poisson","data","motif","wrc", and "all". Specifying
#' "poisson" will result in mutations that can occur anywhere in the heavy chain region, with each nucleotide having an equal probability
#' for a mutation event. Specifying "data" focuses mutation events during SHM in the CDR regions (based on IMGT), and
#' there will be an increased probability for transitions (and decreased probability for transversions). Specifying
#' "motif" will cause neighbor dependent mutations based on a mutational matrix from high throughput sequencing
#' data sets (Yaari et al., Frontiers in Immunology, 2013). "wrc" allows for only the WRC mutational hotspots
#' to be included (where W equals A or T and R equals A or G). Specifying "all" will use all four types of mutations
#' during SHM branching events, where the weights for each can be specified in the "SHM.nuc.prob" parameter.
#' @param SHM.nuc.prob Specifies the rate at which nucleotides change during speciation (SHM) events. This parameter
#' depends on the type of mutation specified by SHM.method. For both "poisson" and "data", the input value determines the probability
#' for each site to mutate (the whole sequence for "poisson" and the CDRs for "data"). For either "motif" or "wrc", the number of
#' mutations per speciation event should be specified. Note that these are not probabilities, but the number of mutations
#' that can occur (if the mutation is present in the sequence). If "all" is specified, the input should be a vector
#' where the first element controls the poisson style mutations, second controls the "data", third controls the "motif"
#' and fourth controls the "wrc".
#' @param baseline.mut Specifies the probability (gamma) for each nucleotide to be mutated inbetween speciation events. These
#' mutations do not cause any branching events. This parameter gives each site a probability to be mutated
#' (in all current sequences) at each time step. Currently these are only Poisson distributed but future
#' releases will change it to allow for other mutation methods.
#' @param SHM.branch.param Describes the probability of undergoing SHM events. This parameter is responsible for
#' describing how likely each sequence will undergo branching events in the phylogeny. The following options are
#' possible: "identical", "uniform", "exponential" ("exp"), "lognormal" ("lognorm"), "normal" ("norm").
#' @param SHM.branch.prob Specifies the probability for a given sequence to undergo SHM events (thus, branching events)
#' This parameter corresponds to the distribution specified in "SHM.branch.prob". For "identical" only one value
#' should be supplied. For "uniform", a vector of length 3 should be specified corresponding to n,min,max respectively
#' (stats::runif(n, min = 0, max = 1)). For "exponential", a single value controlling the rate parameter (from stats::rexp()) should be supplied. For "lognorm"
#' a vector of length two should be supplied, with the first value corresponding to meanlog and the second
#' corresponding to sdlog (from stats::rlnorm). Similarly, for "normal" distribution, two values corresponding to
#' the mean and standard deviation (respectively) should be supplied.
#' @param species Either "mus" for C57BL/6 germline genes or "hum" for human germline genes. These genes were
#' taking from IMGT. When more than one allele was present for a given gene, the first was used.
#' @param max.VDJ The maximum number of VDJ events allowed. These VDJ events are independent of each other
#' but use the same VDJ segments to create a new branching event in the tree at the unmutated germline.
#' @param VDJ.branch.prob The probabilty of a new VDJ recombination event of occuring. For the singleLineage function
#' this will result in a branching event at the site of the unmutated germline. For fullRepertoire function, this
#' will cause a new tree to begin.
#' @param proportion.sampled Value ranging from 0 and 1 specifying the proportion of sequences to be sampled at each time point.
#' Specifiying 1 indicates that all sequences will be recovered at each time point, whereas 0.5 will sample half of the
#' sequences.
#' @param sample.time Integer array indicating the time points at which sampling events should occur.
#' @param chain.type String determining whether heavy or light chain should
#' be simulated. Either "heavy" for heavy chains or "light" for light
#' chains. Heavy chains will have V-D-J recombination, whereas light chain
#' will just have V-J recombination.
#' @param vdj.model Specifies the model used to simulate V-D-J recombination. Can be
#' either "naive" or "data". "naive" is chain independent and does not differentiate
#' between different species. To rely on the default "experimental" options, this
#' should be "data" and the parameter vdj.insertion.mean should be "default". This will
#' allow for different mean additions for either the VD and JD junctions and will
#' differ depending on species.
#' @param vdj.insertion.mean Integer value describing the mean number of nucleotides to be inserted during
#' simulated V-D-J recombination events. If "default" is entered, the mean will be normally distribut
#' @param vdj.insertion.stdv Integer value describing the standard deviation corresponding to insertions
#' of V-D-J recombination. No "default" parameter currently supported but will be updated
#' with future experimental data. This should be a number if using a custom distribution
#' for V-D-J recombination events, but can be "default" if using the "naive" vdj.model
#' or the "data", with vdj.insertion.mean set to "default".
#' @return Returns a nested list containing both sequence information and phylogenetic trees.
#' If "output" is the returned object, then output[[1]][[1]] is an array of the simulated sequences
#' output[[2]][[1]] is an array names corresponding to each sequence. For example, output[[2]][[1]][1]
#' is the name of the sequence corresponding to output[[1]][[1]][1]. The simulated tree of this is found in
#' output[[3]][[1]]. The length of the output list is determined by the number of sampling points
#' Thus if you have two sampling points, output[[4]][[1]] would be a character array holding the sequences
#' with output[[5]][[1]] as a character array holding the corresponding names. Then the sequences recovered
#' second sampling point would be stored at output[[6]][[1]], with the names at output[[7]][[1]]. This
#' nested list was designed for full antibody repertoire simulations, and thus, may seem unintuitive
#' for the single lineage function. The first sequence and name corresponds to the germline sequence
#' that served as the root of the tree.
#' @seealso fullRepertoire
#' @export
#' @examples
# singleLineage(max.seq.num=40,max.timer=150,
#  SHM.method="naive",SHM.nuc.prob = 15/350,
#  baseline.mut = 0.0008, SHM.branch.prob = "identical",
#  SHM.branch.param = 0.05, species="mus",
#  max.VDJ = 1, VDJ.branch.prob = 0.1,
#  proportion.sampled = 1, sample.time = 50,
#  chain.type="heavy",vdj.model="naive",
#  vdj.insertion.mean=4,vdj.insertion.stdv=2)

singleLineage <- function(max.seq.num,
                          max.timer,
                          SHM.method,
                          SHM.nuc.prob,
                          baseline.mut,
                          SHM.branch.prob,
                          SHM.branch.param,
                          species,
                          max.VDJ,
                          VDJ.branch.prob,
                          proportion.sampled,
                          sample.time,
                          chain.type,
                          vdj.model,
                          vdj.insertion.mean,
                          vdj.insertion.stdv){
  max.seq.num <- max.seq.num + 1
  output_list <- list()
  max_SHM <- max.seq.num - max.VDJ
  if(chain.type=="heavy" && species=="mus" || species=="mouse" || species=="blc6"){
    indV <- sample(x = 1:nrow(ighv_mus_df),size = 1,replace = FALSE)
    indD <- sample(x = 1:nrow(ighd_mus_df),size=1,replace=FALSE)
    indJ <- sample(x = 1:nrow(ighj_mus_df),size=1,replace=FALSE)
    vseq <- as.character(ighv_mus_df$seq[indV])
    dseq <- as.character(ighd_mus_df$seq[indD])
    jseq <- as.character(ighj_mus_df$seq[indJ])
    germline <- paste(as.character(ighv_mus_df[[2]][indV]),as.character(ighd_mus_df[[2]][indD]),as.character(ighj_mus_df[[2]][indJ]),sep="")
    germline_v <- as.character(ighv_mus_df[[1]][indV])
    germline_d <- as.character(ighd_mus_df[[1]][indD])
    germline_j <- as.character(ighj_mus_df[[1]][indJ])
  }
  else if(chain.type=="heavy" && species=="hum" || species=="human"){
    indV <- sample(x = 1:nrow(ighv_hum_df),size = 1,replace = FALSE)
    indD <- sample(x = 1:nrow(ighd_hum_df),size=1,replace=FALSE)
    indJ <- sample(x = 1:nrow(ighj_hum_df),size=1,replace=FALSE)
    vseq <- as.character(ighv_hum_df[[2]][indV])
    dseq <- as.character(ighd_hum_df[[2]][indD])
    jseq <- as.character(ighj_hum_df[[2]][indJ])
    germline <- paste(as.character(ighv_hum_df[[2]][indV]),as.character(ighd_hum_df[[2]][indD]),as.character(ighj_hum_df[[2]][indJ]),sep="")
    germline_v <- as.character(ighv_hum_df[[1]][indV])
    germline_d <- as.character(ighd_hum_df[[1]][indD])
    germline_j <- as.character(ighj_hum_df[[1]][indJ])
  }
  else if(chain.type=="light" && species=="mus" || species=="mouse" || species=="blc6"){
    indV <- sample(x = 1:nrow(rbind(iglv_mus_df,igkv_mus_df)),size = 1,replace = FALSE)
    #indD <- sample(x = 1:nrow(rbind(igld_mus_df,igkd_mus_df)),size=1,replace=FALSE)
    indJ <- sample(x = 1:nrow(rbind(iglj_mus_df,igkj_mus_df)),size=1,replace=FALSE)
    vseq <- as.character(rbind(iglv_mus_df,igkv_mus_df)$seq[indV])
    dseq <- ""
    jseq <- as.character(rbind(iglj_mus_df,igkj_mus_df)$seq[indJ])
    germline <- paste(as.character(rbind(iglv_mus_df,igkv_mus_df)[[2]][indV]),as.character(rbind(iglj_mus_df,igkj_mus_df)[[2]][indJ]),sep="")
    germline_v <- as.character(rbind(iglv_mus_df,igkv_mus_df)[[1]][indV])
    germline_d <- ""
    germline_j <- as.character(rbind(iglj_mus_df,igkj_mus_df)[[1]][indJ])
  }
  else if(chain.type=="light" && species=="hum" || species=="human"){
    indV <- sample(x = 1:nrow(rbind(iglv_hum_df,igkv_hum_df)),size = 1,replace = FALSE)
    #indD <- sample(x = 1:nrow(rbind(igld_hum_df,igkd_hum_df)),size=1,replace=FALSE)
    indJ <- sample(x = 1:nrow(rbind(iglj_hum_df,igkj_hum_df)),size=1,replace=FALSE)
    vseq <- as.character(rbind(iglv_hum_df,igkv_hum_df)$seq[indV])
    dseq <- ""
    jseq <- as.character(rbind(iglj_hum_df,igkj_hum_df)$seq[indJ])
    germline <- paste(as.character(rbind(iglv_hum_df,igkv_hum_df)[[2]][indV]),as.character(rbind(iglj_mus_df,igkj_mus_df)[[2]][indJ]),sep="")
    germline_v <- as.character(rbind(iglv_hum_df,igkv_hum_df)[[1]][indV])
    germline_d <- ""
    germline_j <- as.character(rbind(iglj_hum_df,igkj_hum_df)[[1]][indJ])
  }


  #germline <- paste(as.character(mus_hum_df[[1]][indV]),as.character(mus_hum_df[[1]][indD]),as.character(mus_hum_df[[1]][indJ]),sep="")
  #germline <- paste(vseq,dseq,jseq,sep="")
#  germline_v <- as.character(mus_hum_df$gene[mus_hum_df$seq==vseq])
#  germline_d <- as.character(mus_hum_df$gene[mus_hum_df$seq==dseq])
#  germline_j <- as.character(mus_hum_df$gene[mus_hum_df$seq==jseq])
  current_seq_count <- 0
  VDJ_count <- 0
  SHM_count <- 0
  next_node <- 1
  output_list <- list()
  if(class(sample.time)=="numeric"){
    for(i in 1:length(sample.time)+3){
      output_list[[i]] <- list()
    }
  }
  seq_array <- c()
  seq_names <- c()
  seq_type <- c()
  germline_name <- paste(germline_v[1],germline_d[1],germline_j[1],sep="")
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
  for(i in 1:max.timer){
    if(length(seq_array)>=max.seq.num) break
    if(i==sample.time[sample_index] && length(seq_array) >= 1){
      current_len <- length(seq_array)
      holding_size <- as.integer(proportion.sampled * current_len)
      holding_ind <- sample.int(n = current_len,
                                size=holding_size,
                                replace = FALSE,
                                prob = rep(x = 1/current_len,current_len))
      output_list[[2*sample_index+2]] <- seq_array[holding_ind]
      output_list[[2*sample_index+3]] <- paste(seq_names[holding_ind],sample.time[sample_index],sep="_")
      if(length(sample.time) > sample_index) sample_index <- sample_index + 1
    }
    if(length(seq_array)>=2){
      seq_array <- .applyBaseLine(seq_array,baseline.mut) #change to skip1
    }
    is_new_VDJ <- sample(x=c(0,1), replace=TRUE, prob=c(VDJ.branch.prob,
                                                        1- VDJ.branch.prob))
    if(is_new_VDJ==0 && current_seq_count<max.seq.num && VDJ_count<max.VDJ){
      seq_array[next_node] <- .VDJ_RECOMBIN_FUNCTION(vseq, dseq,jseq,
                                                     method=vdj.model,
                                                     chain.type=chain.type,
                                                     species=species,
                                                     vdj.insertion.mean=vdj.insertion.mean,
                                                     vdj.insertion.stdv=vdj.insertion.stdv)
      seq_names[next_node] <- paste(paste("L",next_node, sep=""),i,sep="_")
      if(next_node>=2){
        output_tree_text <- .branchingProcess3(output_tree_text,0,next_node,"lineage")
      }

      current_seq_count <- current_seq_count + 1
      VDJ_count <- VDJ_count + 1
      next_node <- next_node + 1
    }
    if(next_node>= max.seq.num) break
    if(current_seq_count>=1){
      for(j in 1:current_seq_count){
        is_new_SHM <- sample(x=c(0,1),1, replace=TRUE, prob=c(new_SHM_prob[next_node], 1- new_SHM_prob[next_node]))
        if (is_new_SHM==0 && next_node<max.seq.num && SHM_count<max_SHM ){
          seq_array[next_node] <- .SHM_FUNCTION_SEQUENCE4(seq_array[j],
                                                         SHM.method,vseq,dseq,jseq,SHM.nuc.prob)
          seq_names[next_node] <- paste(paste("S",next_node, sep=""),i,sep="_")
          output_tree_text <- .branchingProcess3(output_tree_text,j,next_node,"SHM")
          next_node <- next_node + 1
          current_seq_count <- current_seq_count + 1
          SHM_count <- SHM_count + 1

        }
        if(length(seq_array)>=max.seq.num) break
      }
    }
    if(length(seq_array)>=max.seq.num) break

  }
  output_list[[1]][[1]] <- append(germline, seq_array)
  output_list[[2]][[1]] <- append(germline_name, seq_names)
  temp.tree <- ape::read.tree(text=output_tree_text)
  temp.tree$tip.label[1] <- germline_name[1]
  output_list[[3]][[1]] <- temp.tree
  for(i2 in 2:length(output_list[[3]][[1]]$tip.label)){
    output_list[[3]][[1]]$tip.label[1] <- output_list[[2]][[1]][1]
    if(length(output_list[[2]][[1]][grep(output_list[[2]][[1]],pattern=paste("S",output_list[[3]][[1]]$tip.label[i2],"_",sep=""))])>0){
      output_list[[3]][[1]]$tip.label[i2] <- output_list[[2]][[1]][grep(output_list[[2]][[1]],pattern=paste("S",output_list[[3]][[1]]$tip.label[i2],"_",sep=""))]
    }
    else if(length(output_list[[2]][[1]][grep(output_list[[2]][[1]],pattern=paste("L",output_list[[3]][[1]]$tip.label[i2],"_",sep=""))])>0){
      output_list[[3]][[1]]$tip.label[i2] <- output_list[[2]][[1]][grep(output_list[[2]][[1]],pattern=paste("L",output_list[[3]][[1]]$tip.label[i2],"_",sep=""))]
    }
  }
  return(output_list)
}



