#!/usr/bin/env Rscript

# Read in quietly to avoid outputting "Loading required package: Rcpp" to stderr.
library(castor, quietly = TRUE)

# Read in ape as well for "read.tree" function.
library(ape)

Args <- commandArgs(TRUE)

# Read in command-line arguments.
full_tree <- read.tree(Args[1])
trait_values <- read.delim(Args[2], check.names=FALSE, row.names=1)
hsp_method <- Args[3]
calc_nsti <- as.logical(Args[4])
calc_ci <- as.logical(Args[5])
check_input_set <- as.logical(Args[6])
weight_setting <- as.logical(Args[7])
num_threads <- Args[8]
predict_outfile <- Args[9]
ci_outfile <- Args[10]

# Function to get CIs for certain HSP methods.
ci_95_states2values <- function(state_probs, states2values, number_of_tips) {
  state_prob_cumsum <- t(apply(state_probs[1:number_of_tips,], 1, cumsum))
  ci_5 <- apply(state_prob_cumsum, 1, function(x) { states2values[min(which(x >= 0.05)),] } )
  ci_95 <- apply(state_prob_cumsum, 1, function(x) { states2values[min(which(x >= 0.95)),] } )
  return(paste(as.character(ci_5), as.character(ci_95), sep = "|"))
}

# If tips in tree contain '_' then read.tree() places single quotes around these tip labels.
# This then causes sorting errors below since the rownames are different between the trait table and the tree.
# Fix this by putting quotes around any labels in the trait table that have a '_'.
for(i in grep("_",rownames(trait_values))) {
 rownames(trait_values)[i] <- paste("'", rownames(trait_values)[i], "'", sep="")
}

# Order the trait table to match the tree tip labels. Set all tips without a value to be NA.
unknown_tips <- full_tree$tip.label[which(! full_tree$tip.label %in% rownames(trait_values))]

unknown_df <- as.data.frame(matrix(NA,
                                  nrow=length(unknown_tips),
                                  ncol=ncol(trait_values)))

rownames(unknown_df) = unknown_tips
colnames(unknown_df) = colnames(trait_values)

trait_values_ordered <- rbind(trait_values, unknown_df)
trait_values_ordered <- trait_values_ordered[full_tree$tip.label, , drop=FALSE]

num_tip <- nrow(trait_values_ordered)

tmp <- find_nearest_tips(full_tree, 
                         target_tips=tip_range[-2],
                         check_input=check_input_set)

if (hsp_method == "pic") {
  
  predict_out <- apply(trait_values_ordered,
                       2,
                       hsp_independent_contrasts,
                       tree=full_tree,
                       weighted=weight_setting,
                       check_input=check_input_set)
  
  predicted_values <- sapply(predict_out[num_tip,], function(x) { x$states })

} else if (hsp_method == "sqp") {

  predict_out <- apply(trait_values_ordered,
                       2,
                       hsp_squared_change_parsimony,
                       tree=full_tree,
                       weighted=weight_setting,
                       check_input=check_input_set)
  
  predicted_values <- sapply(predict_out[num_tip,], function(x) { x$states })
  
} else if(hsp_method == "mk_model" | hsp_method == "emp_prob" | hsp_method = "mp") {

  # Typically to convert from state index to value you would just need to substract 1.
  # However, a user could also use this script without integer data starting at 0 so, the below
  # commands will allow for that.
  trait_states_mapped_full <- apply(trait_values_ordered,
                                    2,
                                    map_to_state_space,
                                    fill_gaps = TRUE,
                                    sort_order = "natural",
                                    include_state_values = TRUE)
  
  trait_states_mapped <- lapply(trait_states_mapped_full, function(x) { x$mapped_states })
  trait_states_values <- lapply(trait_states_mapped_full, function(x) { x$state_values })
  
  if(hsp_method == "mk_model") {
    hsp_out_models <- lapply(trait_states_mapped,
                            hsp_mk_model,
                            tree = full_tree,
                            rate_model = "SRD",
                            check_input = check_input_set,
                            Nthreads = num_threads)
                            
  } else if (hsp_method == "emp_prob") {
    
    hsp_out_models <- lapply(trait_states_mapped,
                            hsp_empirical_probabilities,
                            tree = full_tree,
                            check_input = check_input_set)
    
  } else if (hsp_method == "mp") {

    hsp_out_models <- lapply(trait_states_mapped,
                            hsp_max_parsimony,
                            tree = full_tree,
                            check_input = check_input_set,
                            transition_costs = "sequential",
                            weight_by_scenarios = TRUE)
                                       
  }
  
  # Get state with highest probability in each case.
  predicted_states <- sapply(hsp_out_models,
                             function(x) { max.col(x$likelihoods[1:num_tip,]) })
  
  state2val <- data.frame(row.names=1:trait_states$Nstates,
                          value=trait_states$state_values)
  
  predicted_values <- apply(predicted_states, 2, function(x) { state2val[x,] } )

  # If calc_ci set then figure out what the assigned trait would be at the 95% CI and output resulting matrix.
  if(calc_ci) {
    ci_values <- data.frame(sapply(hsp_out_models, function(x) { ci_95_states2values(x$likelihoods, state2val, num_tip) }))
    ci_values$tips <- full_tree$tip.label
    ci_values <- ci_values[, c("tips", colnames(trait_values_ordered))]
    write.table(ci_values, file=ci_outfile, sep="\t", quote=FALSE, row.names=FALSE)
   }
      
}

# Add tips as first column of predicted_values.
predicted_values <- data.frame(predicted_values)
predicted_values$tips <- full_tree$tip.label
predicted_values <- predicted_values[, c("tips", colnames(trait_values_ordered))]

# Calculate NSTI per tip and add to output as last column if option set.
if(calc_nsti) {
  tip_range <- 1:num_tip
  predicted_values$nsti <- sapply(tip_range, 
                                  function(x) { find_nearest_tips(full_tree, 
                                                                  target_tips=tip_range[-x],
                                                                  check_input=check_input_set)$nearest_distance_per_tip[x]})
}

# Write out predicted values.
write.table(predicted_values, file=predict_outfile, row.names=FALSE, quote=FALSE, sep="\t")
