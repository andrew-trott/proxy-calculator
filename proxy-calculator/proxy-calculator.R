# THE IDEA BEHIND THIS CODE IS TO STREAMLINE PROXY CALCULATIONS, LIMITING THE NUMBER OF EDITS NEEDED FOR THE CODE TO RUN.
# THE PURPOSE BEHIND EACH LINE OF CODE IS EXPLAINED AS MUCH AS POSSIBLE.
# Once options have been set and data has been uploaded in relative abundance form, run the script in its entirety, the user options will take care of the rest.



#Set your working directory (you might need to change any \ to /)
#setwd("C:/Users/etc...") #This is the working directory you wish to save the final data to.


#USER OPTIONS
data  <- CH_ES_DATA ###########################Just put the name of your uploaded data set here.
use_parallel <- TRUE ###############################TRUE for parallel, FALSE for sequential. Parallel execution distributes the workload to all processors on your system, allowing the code to run far faster. Running parallel will mean that you won't see a progress bar.
use_log_ratio <- FALSE #############################TRUE for log(numerator) - log(denominator), FALSE for numerator/denomiator. This only changes the way in which the ratio is calculated.
max_k <- 5 #########################################This number limits the number of combinations allowed for numerator and denominator (this can significantly speed up the calculation process)
variable_cols <- data[,2:21] #######################These are the columns the proxy variables occupy in your data set (data[,column:to this column])
environmental_col <- log(data$`ca mg/L`) ###########This is the environmental parameter you wish to test against (data$["insert parameter here"]). Put log(data$parameter) if you want to test against the log of the parameter.
env_filter_min   <- NA       #######################This is the lower limit of the environmental parameter. Leave it as NA if you don't want to apply one. If you only want either an upper or lower bound, just set one and leave the other as NA.
env_filter_max   <- NA      ########################This is the upper limit of the environmental parameter. Again, leave empty if you don't want to apply one. 
r_threshold <- 0.60  ###############################This is the threshold for the r-values that you are interested in, change this as necessary.


#CHECKS FOR USERS
head(variable_cols) #These are all for you to double check that the data looks good, and no problems have crept in.
head(environmental_col)

#CHECKS FOR THE SCRIPT
stopifnot(is.data.frame(data)) #This section is the programme checking to make sure the data looks suitable and, if it's not, it'll stop early to avoid you wasting time.
stopifnot(ncol(variable_cols) > 0, is.numeric(environmental_col))
if (any(is.na(variable_cols))) warning("Missing values found in proxy variables") #These lines won't stop the code, but they will warn you if there are issues, so keep an eye on the console.
if (any(is.na(environmental_col))) warning("Missing values found in environmental parameter")

#BLOCK 1: Set the limits
if (!is.na(env_filter_min) || !is.na(env_filter_max)) { #This checks to see if you set any limits to the environmental parameter. If you didn't then this whole thing is skipped.
  keep <- rep(TRUE, length(environmental_col)) #This makes a logical vector full of Trues that the later parts are applied to.
  if (!is.na(env_filter_min)) keep <- keep & (environmental_col >= env_filter_min) #Assuming you set a lower boundary, this is applied here.
  if (!is.na(env_filter_max)) keep <- keep & (environmental_col <= env_filter_max) #This is the same just for the upper boundary.
  
  variable_cols     <- variable_cols[keep, , drop = FALSE] #This creates a subset of the paramters, including only those that passed the boundary checks.
  environmental_col <- environmental_col[keep] #This then applies that to the environmental data.
  
  message("Kept ", sum(keep), " rows where ", #This section sends a message in the console to let you know how many rows were kept and to reiterate the boundaries.
          ifelse(is.na(env_filter_min), "-∞", env_filter_min),
          " ≤ env_col ≤ ",
          ifelse(is.na(env_filter_max), "+∞", env_filter_max))
}


# BLOCK 2: Check for packages, and install and load any that are missing.
use_parallel <- isTRUE(use_parallel) #this runs this section of code if you selected the option to use parallel execution. This also defaults the use_parallel setting to FALSE if it is anything other than TRUE.

pkgs_needed <- c("vegan", "sampling", "progress", "doParallel") #this sets up the packages that are needed.
if (use_parallel) { 
  pkgs_needed <- c(pkgs_needed, "doParallel", "foreach") #this grabs the packages that are needed for parallel execution and adds them to the list of required packages.
}

pkgs_to_install <- setdiff(pkgs_needed, rownames(installed.packages())) #this defines the difference between the packages that you have already installed, and those you still need. This saves you from reinstalling packages you already have.
if (length(pkgs_to_install) > 0) { #assuming you need to install packages, this will run.
  install.packages(pkgs_to_install, dependencies = TRUE) #grabs and installs the required packages.
}

suppressPackageStartupMessages({ #just stops the messages from showing up in the console.
  lapply(pkgs_needed, function(pkg) library(pkg, character.only = TRUE)) #this loads the required packages *quietly*
})

#BLOCK 3
#ParallelExecution prep. This will allow the code to run in parallel, allowing for multiple CPU cores to work simultaneously on a problem. When run, there may be warning messages to do with `closing unused connections`, these are okay and don't impact the code.
if (use_parallel) {
  num_cores <- detectCores() - 1  # Leave one free! #PE
  cl <- makeCluster(num_cores) #creates a cluster using the cores defined above. #PE
  registerDoParallel(cl) #Registers the defined cluster for parallel execution. #PE
  on.exit(stopCluster(cl)) #this ensures that, if the code quits it is completed, the cluster is stopped #if you stop it yourself, run the stopCluster command manually!
  message("Using ", num_cores, " cores :D") #this will allow you to see, as the code is running, how many cores it is using. 
} else {
  message("Parallel mode: OFF. Running sequentially.")
}


#BLOCK 4
#convert to matrix, run the 0 replacer, and convert back again.
my.data.values <-variable_cols #rename for clarity, not strictly necessary
my.mat <- as.matrix(my.data.values) #converts the dataframe to a matrix, to allow for 0 replacement. If it were a dataframe, this step may fail (seems to want to change entire rows that have 0s)
zero.count <- sum(my.mat == 0) #find these 0s
my.mat[my.mat == 0] <- runif(zero.count, 0, 0.0002) #set them to a very small, but still >0, value
my.data.values <- as.data.frame(my.mat) #convert back to a data frame for the next steps.


#BLOCK 5
#preparation for proxy calculation code
combinationsN <- function(Data, max_k) { #this grabs the max_k value from earlier
  n <- ncol(Data) #sets n to be equal to the column number of your data set.
  stopifnot(max_k <= n) #stops the code if there has been an error, and the data does not have the expected parameters.
  Combinations_list <- list() #creates a list to show the included variables in binary form (this means that an included variable is shown as 1, the rest will be 0)
  Sums_list <- list() #creates an empty list to store the row-wise sum of each variable.
  j <- 0 #this tracks how many combinations have been seen and uses it to store the data.
  for (k in 1:max_k) { #this loops through all combination of 1 variable, through to max_k variables (e.g. 6)
    combos_k <- combn(n, k, simplify = FALSE) #this generates all the combinations from 1 to max_k variables.
    for (combo in combos_k) { #for all of these combinations...
      j <- j + 1 #...first it adds a tick to the counter...
      Combinations_list[[j]] <- as.integer(seq_len(n) %in% combo) #...then it creates the binary vector for the list...
      sum_vals <- if (length(combo) == 1) { #...then it sums those vectors, if there is more than one variable.
        Data[, combo]
      } else {
        rowSums(Data[, combo, drop = FALSE])
      }
      Sums_list[[j]] <- sum_vals #this then saves the list.
    }
  }
  Combinations <- do.call(rbind, Combinations_list)  #this stacks all the binary vectors, row-by-row, into a matrix (#combos × n)
  Sums         <- do.call(cbind, Sums_list)          # (samples × #combos)
  list(Sums = Sums, Combinations = Combinations)     # this returns a named list with the two matrices.
}

#for n samples, there are 2^n -1 combinations, this should appear where combination no. appears below.
res_comb <- combinationsN(my.data.values, max_k) #generates all the combinations up to max_k variables. Required later for the calculations.
if (any(res_comb$Sums < 1e-6)) warning("Some combinations have very small sums; this may lead to unstable ratio values.") #This line just sends you a warning message in the unlikely situation that some of the denominators are 0.
n_combos  <- ncol(res_comb$Sums) #grabs the number of combinations.
cat("Generated", n_combos, "combinations (each up to", max_k, "variables).\n") #these are message prompts for your convenience!
cat("Dims of Sums:", dim(res_comb$Sums), "\n") #These two lines are again just messages for you, to let you know how many combinations were generated and the dimensions of the matrix.
cat("Dims of Combinations:", dim(res_comb$Combinations), "\n")


#BLOCK 6
#this is the calculation block, what happens here depends on whether or not you've used parallel execution and on the ratio type that you want.
if (isTRUE(use_parallel)) { #again checks your response to the parallel execution setting.
  
  # Parallel version.
  R10_list <- foreach(j = 1:n_combos, #this section runs the loop in parallel using %dopar%
                      .combine = rbind, 
                      .packages = "stats") %dopar% {
                        pR  <- numeric(n_combos) #these lines create a numeric vector to store the correlation values in order to speed up calculations.
                        num <- pmax(res_comb$Sums[, j], 1e-6) #this takes the numerator values and replaces any very small or 0 values with 1e-6. This is just to avoid errors.
                        
                        for (i in 1:n_combos) { #this loops through all denominators.
                          if (i == j) {  #this line skips any ratio that is simply A/A. These will always be 1, of course, and though it doesn't do much (like 0.1% faster) but it will still speed things up a bit!
                            message(paste("Skipping self-comparison at combo", i))
                            next
                          }
                          den <- pmax(res_comb$Sums[, i], 1e-6) #same process that was used for the numerators.
                          
                          if (isTRUE(use_log_ratio)) { #This is the log ratio section, and is used depending on your choice from the settings.
                            # log(numerator) - log(denominator) 
                            pR[i] <- suppressWarnings(cor(environmental_col, log(num) - log(den))) #calculates the correlation between the environmental parameter and the log - log ratios.
                          } else { #else here means that you put FALSE in the log ratio setting, therefore running the n/d ratio instead.
                            # numerator / denominator
                            pR[i] <- suppressWarnings(cor(environmental_col, num / den)) #calculates the correlation between the environmental parameter and the n/d ratios.
                          }
                        }
                        
                        ord <- order(pR, decreasing = TRUE, na.last = TRUE) #this orders the correlations from highest to lowest.
                        if (length(ord) >= 10) { #this section picks the top 10 values, or pads the matrix with NAs if there are less than 10. 
                          ord[1:10] #this gives the best 10 values from the "ord" list
                        } else {
                          c(ord, rep(NA_integer_, 10 - length(ord))) #this makes sure the list has at least 10 entries, just in case something happens to stop that.
                        }
                      }
  
  R10 <- R10_list #this just assigns the result to an output variable.
  
} else { #this is the same as the previous version, except it is run in sequence, not parallel.
  
  # Sequential version with a progress bar for your convinience.
  R10 <- matrix(NA_integer_, nrow = n_combos, ncol = 10) #this creates a matrix with the correct dimensions.
  pb <- progress_bar$new( #this is the section that creates the progress bar!
    format = "  [:bar] :percent ETA: :eta", 
    total  = n_combos, clear = FALSE, width = 60
  )
  
  for (j in 1:n_combos) { #This section...
    pR  <- numeric(n_combos)
    num <- pmax(res_comb$Sums[, j], 1e-6)
    
    for (i in 1:n_combos) {
      if (i == j) {   
        message(paste("Skipping self-comparison at combo", i))
        next
      }
      den <- pmax(res_comb$Sums[, i], 1e-6)
      
      if (isTRUE(use_log_ratio)) {
        pR[i] <- suppressWarnings(cor(environmental_col, log(num) - log(den)))
      } else {
        pR[i] <- suppressWarnings(cor(environmental_col, num / den))
      }
    }
    
    ord <- order(pR, decreasing = TRUE, na.last = TRUE)
    if (length(ord) >= 10) {
      R10[j, ] <- ord[1:10]
    } else {
      R10[j, ] <- c(ord, rep(NA_integer_, 10 - length(ord))) #...does the same as the parallel one.
    }
    
    pb$tick() #this is part of the progress bar, each time the loop reaches here it ticks it.
  }
}

#BLOCK 7
#this section is the export prep section, and will take a while too unfortunately.
var_names <- colnames(my.data.values) #this gets the original names from the data set that you uploaded.
top_results_raw <- NULL #this makes an object frame to hold the results before we make the final data frame.
if (isTRUE(use_parallel)) { #checking your user options.
  top_results_raw <- foreach(j = 1:n_combos, .combine = rbind, .packages = "stats") %dopar% { #uses the parallel execution to loop over each numerator.
    if (!use_log_ratio) num <- pmax(res_comb$Sums[, j], 1e-6) #again limits the lowest values to avoid log(0) errors.
    
    results_list <- list() #this makes an empty list to collect one data frame per significant combination
    
    for (i in R10[j, ]) { #loops the top 10 denominator indices for the numerators
      if (is.na(i)) next #this just makes it skip any NAs
      
      if (use_log_ratio) { #This section.....
        log_j <- log(pmax(res_comb$Sums[, j], 1e-6))
        log_i <- log(pmax(res_comb$Sums[, i], 1e-6))
        ratio <- log_j - log_i
      } else {
        den   <- pmax(res_comb$Sums[, i], 1e-6)
        ratio <- num / den
      } #.....uses either the log or division ratio, your choice, and recomputes their correlation to be stored in the data frame.
      
      r_val <- suppressWarnings(cor(environmental_col, ratio)) #this figures out the Pearson correlation of the variables, and silences warnings.
      if (is.na(r_val) || abs(r_val) < r_threshold) next #if the results are NAs, this just skips them and so speeds up the next steps.
      
      df     <- length(environmental_col) - 2 #this is the degrees of freedom for the correlation test.
      t_val  <- r_val * sqrt(df / (1 - r_val^2)) #this transforms the Pearson correlation (r) into a t-statistic (using the equation given) with the idea being that a large t-value means a significant correlation. We need this to compute the p-value.
      p_val  <- 2 * pt(-abs(t_val), df = df) #this is the equation to calculate the p-value
      r2_val <- r_val^2 #this gives the r^2 value
      
      combo_j <- res_comb$Combinations[j, ] #This section...
      combo_i <- res_comb$Combinations[i, ]
      
      numer_names <- var_names[combo_j == 1]
      denom_names <- var_names[combo_i == 1] #...extracts the indicator rows and converts them to readable names for the formula column.
      
      formula_str <- if (use_log_ratio) { #checks your user options and builds a representation of the formula using that choice.
        paste0("log(", paste(numer_names, collapse = " + "), ") - log(", paste(denom_names, collapse = " + "), ")") #for log - log.
      } else {
        paste0("(", paste(numer_names, collapse = " + "), ") / (", paste(denom_names, collapse = " + "), ")") #for simple division.
      }
      
      results_list[[length(results_list) + 1]] <- data.frame( #This section...
        Numerator         = j,
        Denominator       = i,
        r                 = r_val,
        r_squared         = r2_val,
        p                 = p_val,
        Numerator_Combo   = paste(combo_j, collapse = ""),
        Denominator_Combo = paste(combo_i, collapse = ""),
        Formula           = formula_str,
        stringsAsFactors  = FALSE #... creates a data frame with all the right names and sets it to the results list.
      )
    }
    
    if (length(results_list) == 0) return(NULL) #if nothing passes the threshold then it makes it NULL so the binder skips it.
    do.call(rbind, results_list) #this binds all single-row data frames in the results list in a single multi-row data frame.
  }
  
} else {
  # Sequential version with progress bar. Otherwise it is exactly the same as the above one.
  top_results_raw <- list()
  pb <- progress_bar$new(
    format = "  [:bar] :percent ETA: :eta",
    total  = n_combos, clear = FALSE, width = 60
  )
  
  for (j in 1:n_combos) {
    if (!use_log_ratio) num <- pmax(res_comb$Sums[, j], 1e-6)
    
    results_list <- list()
    
    for (i in R10[j, ]) {
      if (is.na(i)) next
      
      if (use_log_ratio) {
        log_j <- log(pmax(res_comb$Sums[, j], 1e-6))
        log_i <- log(pmax(res_comb$Sums[, i], 1e-6))
        ratio <- log_j - log_i
      } else {
        den   <- pmax(res_comb$Sums[, i], 1e-6)
        ratio <- num / den
      }
      
      r_val <- suppressWarnings(cor(environmental_col, ratio))
      if (is.na(r_val) || abs(r_val) < r_threshold) next
      
      df     <- length(environmental_col) - 2
      t_val  <- r_val * sqrt(df / (1 - r_val^2))
      p_val  <- 2 * pt(-abs(t_val), df = df)
      r2_val <- r_val^2
      
      combo_j <- res_comb$Combinations[j, ]
      combo_i <- res_comb$Combinations[i, ]
      
      numer_names <- var_names[combo_j == 1]
      denom_names <- var_names[combo_i == 1]
      
      formula_str <- if (use_log_ratio) {
        paste0("log(", paste(numer_names, collapse = " + "), ") - log(", paste(denom_names, collapse = " + "), ")")
      } else {
        paste0("(", paste(numer_names, collapse = " + "), ") / (", paste(denom_names, collapse = " + "), ")")
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Numerator         = j,
        Denominator       = i,
        r                 = r_val,
        r_squared         = r2_val,
        p                 = p_val,
        Numerator_Combo   = paste(combo_j, collapse = ""),
        Denominator_Combo = paste(combo_i, collapse = ""),
        Formula           = formula_str,
        stringsAsFactors  = FALSE
      )
    }
    
    if (length(results_list) > 0) {
      top_results_raw[[length(top_results_raw) + 1]] <- do.call(rbind, results_list)
    }
    
    pb$tick() #ticks the progress bar
  }
  
  top_results_raw <- do.call(rbind, top_results_raw)
}

# Ensure a valid data frame is returned
top_results <- if (!is.null(top_results_raw) && nrow(top_results_raw) > 0) { #assuming that there were valid results, this section wraps the raw result frame into the top_results. If there were no valid results, it just makes an empty data frame.
  top_results_raw
} else {
  data.frame() #the empty data frame mentioned above.
  message("⚠ No significant results found (no correlations met the self-selected r and p thresholds).")
}

#BLOCK 8
#export data in a .csv format and voila!
if (nrow(top_results) > 0) {  #This section...
  top_results <- subset(top_results, p < 0.05)
  top_results <- top_results[order(-abs(top_results$r)), ]
  top20 <- head(top_results, 20) #...checks for the results, filters them to ensure statistically significant p-values, sorts the resulting selection by highest r-value, and takes the top 20.
} else {
  message("⚠ No significant results found (none met r and p thresholds).")
  top20 <- data.frame()  # ensures this is defined, if nothing met the criteria then it creates an empty data frame.
}

write.csv(top20, "Proxy Calculation Results.csv", row.names = FALSE) #writes the top20 data frame to a .csv file that it then exports to whichever directory you selected.
message("▶ Done. Results written to Proxy Calculation Results.csv.") #Just to let you know that it was successful :)
