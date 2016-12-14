#' blc6_v_df
#'
#' C57BL/6 germline v-gene segments. When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 164 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"blc6_v_df"

#' blc6_d_df
#'
#' C57BL/6 germline d-gene segments. When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 27 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"blc6_d_df"

#' blc6_j_df
#'
#' C57BL/6 germline j-gene segments. When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 4 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"blc6_j_df"

#' hum_v_df
#'
#' human germline v gene segments. When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 119 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"hum_v_df"


#' hum_d_df
#'
#' human germline v gene segments. When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 37 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"hum_d_df"


#' hum_j_df
#'
#' human germline v gene segments. When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 6 rows and 10 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"hum_j_df"



#' one_spot_df
#'
#' WRC hotspot mutations taken from Yaari et al., Frontiers in Immunology, 2013.
#' These include only the mutations following the WRC pattern,
#' where W equals A or T and R equals A or G). Custom mutation hotspots can be supplied
#' by modifying this dataframe. Repeating particular hotspot entries allows
#' for the hotspot to mutate more than one time per SHM event.
#'
#' @format A data frame with 32 rows and 6 variables:
#' \describe{
#'   \item{pattern}{Character array where each entry corresponds to a 5 base motif. The
#'   mutation probabilities correspond to the middle nucleotide in each 5mer.}
#'   \item{toA}{The probability for the middle nucleotide in "pattern" to mutate to an adenine}
#'   \item{toC}{The probability for the middle nucleotide in "pattern" to mutate to an cytosine}
#'   \item{toG}{The probability for the middle nucleotide in "pattern" to mutate to an guanine}
#'   \item{toT}{The probability for the middle nucleotide in "pattern" to mutate to an thymine}
#'   \item{Source}{The origin of how this motif was discovered. Either Inferred or Experimental}
#' }
#' @source Yaari et al., Frontiers in Immunology, 2013
"one_spot_df"


#' hotspot_df
#'
#' Hotspot mutations taken from Yaari et al., Frontiers in Immunology, 2013.
#' This contains transition probabilities for all 5mer combinations based
#' on high throughput sequencing data. The transition probabilities are for
#' the middle nucleotide in each 5mer set. This can be customized by changing the genes
#' and sequences. Custom mutation hotspots can be supplied
#' by modifying this dataframe. Repeating particular hotspot entries allows
#' for the hotspot to mutate more than one time per SHM event.
#'
#'  @format A data frame with 32 rows and 6 variables:
#' \describe{
#'   \item{pattern}{Character array where each entry corresponds to a 5 base motif. The
#'   mutation probabilities correspond to the middle nucleotide in each 5mer.}
#'   \item{toA}{The probability for the middle nucleotide in "pattern" to mutate to an adenine}
#'   \item{toC}{The probability for the middle nucleotide in "pattern" to mutate to an cytosine}
#'   \item{toG}{The probability for the middle nucleotide in "pattern" to mutate to an guanine}
#'   \item{toT}{The probability for the middle nucleotide in "pattern" to mutate to an thymine}
#'   \item{Source}{The origin of how this motif was discovered. Either Inferred or Experimental}
#' }
#' @source Yaari et al., Frontiers in Immunology, 2013
"hotspot_df"

