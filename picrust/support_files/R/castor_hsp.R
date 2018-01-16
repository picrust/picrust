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
no_round_setting <- as.logical(Args[7])
num_threads <- Args[8]
predict_outfile <- Args[9]
ci_outfile <- Args[10]

# Function to get CIs for certain HSP methods.
ci_95_states2values <- function(state_probs, states2values, number_of_tips) {
  state_prob_cumsum <- t(apply(state_probs[1:number_of_tips, , drop=FALSE], 1, cumsum))
  ci_5 <- apply(state_prob_cumsum, 1, function(x) { states2values[min(which(x >= 0.05)),] } )
  ci_95 <- apply(state_prob_cumsum, 1, function(x) { states2values[min(which(x >= 0.95)),] } )
  return(c(ci_5, ci_95))
}

# Function to identify traits that have only 1 state across tips (for pre-processing mk_model input).
identify_single_state <- function(state_vec) {
  
  state_vec_noNA <- na.exclude(state_vec)
  state_vec_noNA_l <- length(state_vec_noNA)
  
  if(length(which(state_vec_noNA == 1)) == state_vec_noNA_l) {
    return(TRUE)
  } else {
    return(FALSE)
  }
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

if (hsp_method == "pic" | hsp_method == "scp") {
  
  if (hsp_method == "pic") {
    predict_out <- apply(trait_values_ordered,
                         2,
                         hsp_independent_contrasts,
                         tree=full_tree,
                         weighted=TRUE,
                         check_input=check_input_set)
    
  } else if (hsp_method == "scp") {
  
    predict_out <- apply(trait_values_ordered,
                         2,
                         hsp_squared_change_parsimony,
                         tree=full_tree,
                         weighted=TRUE,
                         check_input=check_input_set)
  }

  predicted_values <- sapply(predict_out, function(x) { x$states[1:num_tip] })

  # Round to nearest integer unless no_round set.
  if(! no_round_setting) {
    predicted_values <- round(predicted_values, digits=0)
  }
  
} else if(hsp_method == "mk_model" | hsp_method == "emp_prob" | hsp_method == "mp") {

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
    
    # mk_model throws an error if there is only 1 known state (e.g. if the trait isn't found in any tips).
    # To get around this can just set all predictions to be this known state for these traits.
    single_state_traits <- which(sapply(trait_states_mapped, identify_single_state) == TRUE)
    
    if(length(single_state_traits) > 0) {
      trait_states_mapped_variable <- trait_states_mapped[-single_state_traits]
    } else {
      trait_states_mapped_variable <- trait_states_mapped
    }
    
    if(length(trait_states_mapped_variable) == 0) {
      stop("Stopping - all input traits have the same state value, which won't work with a Markov model.")
    }
    
    hsp_out_models <- lapply(trait_states_mapped_variable,
                             hsp_mk_model,
                             tree = full_tree,
                             rate_model = "SRD",
                             check_input = check_input_set,
                             Nthreads = num_threads)
    
    # If applicable add the invariable traits into the hsp_mk_model output
    if(length(single_state_traits) > 0) {
      for(single_state_trait in names(single_state_traits)) {
        hsp_out_models[[single_state_trait]] <- list(likelihoods = matrix(rep(1, num_tip), nrow=num_tip, ncol=1), success = TRUE)
      }
    }
    
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
  
  max_num_states <- max(sapply(trait_states_mapped_full, function(x) { length(x$state_values)}))
  
  state2val <- sapply(trait_states_mapped_full, function(x) { c(x$state_values, rep(NA, max_num_states-length(x$state_values)))})
  
  predicted_values <- sapply(colnames(predicted_states), function(x) { state2val[predicted_states[,x],x] })

  # If calc_ci set then figure out what the assigned trait would be at the 95% CI and output resulting matrix.
  if(calc_ci) {
    ci_values <- data.frame(sapply(names(hsp_out_models), function(x) { ci_95_states2values(hsp_out_models[[x]]$likelihoods,
                                                                                            state2val[ , x, drop = FALSE],
                                                                                            num_tip) }), check.names = FALSE)
                                                                                                        
    
    ci_values_ci_5 <- ci_values[1:num_tip, , drop=FALSE]
    ci_values_ci_95 <- ci_values[(num_tip+1):(num_tip*2), , drop=FALSE]
    
    colnames(ci_values_ci_5) <- paste(colnames(ci_values_ci_5), "5", sep="_")
    colnames(ci_values_ci_95) <- paste(colnames(ci_values_ci_95), "95", sep="_")
    
    ci_values <- cbind(ci_values_ci_5, ci_values_ci_95)
    
    # Sort column names so that 5% and 95% CIs are next to each other.
    ci_values <- ci_values[ , order(names(ci_values))]
    
    orig_ci_colnames <- colnames(ci_values)
    
    ci_values$tips <- full_tree$tip.label
    ci_values <- ci_values[, c("tips", orig_ci_colnames)]
    write.table(ci_values, file=ci_outfile, sep="\t", quote=FALSE, row.names=FALSE)
   }
}

# Add tips as first column of predicted_values.
predicted_values <- data.frame(predicted_values, check.names = FALSE)
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
