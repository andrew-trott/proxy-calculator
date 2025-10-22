README:

A self‑contained R script to streamline the development and evaluation of proxy ratios.  
Allows you to explore combinations of variables, compute either simple or log‑ratios, and test their correlation against an environmental parameter—all with minimal manual steps.
Once your data has been uploaded, a working directory set, and the user options selected, just run the code in its entirety.

--------------------------------------------

The purpose of this script is to streamline the development of proxy ratios. To make this as flexible as possible, the script has user options that allow you to decide what sort of proxy you want, and to control the variables used for it which will explained below. For ease of use, the code will do everything for you (including installing the appropriate packages) and can be used to calculate any number of variable combinations.

It is worth keeping an eye on the console during the runtime, as there are messages that will pop up to keep you informed about the state of the run and the data. It will typically get to a point and then remain there for the majority of the runtime, but so long as you see that there is a "stop" button in the top righthand corner of the console, the code is still running. If you do not run in parallel, then the code will provide a progress tick bar that will be a more direct way of informing you.

Note: You must have R ≥ 3.6.0 installed. All required packages will be installed automatically

--------------------------------------------
FEATURES
--------------------------------------------

- Flexible ratio types: simple -> ((A + B …)/(C + D …)) or log‑ratio -> (log(A + B …) - log(C + D …)).  
- Automatic package management: installs and loads any missing CRAN packages.  
- Parallel execution: leverages all available CPU cores for huge combinatorial searches. This should be suitable for supercomputer usage.  
- Progress feedback: console messages in parallel mode, and an additional interactive progress bar in sequential mode.  
- Statistical testing: computes Pearson’s (r), t-statistic, p-value, and R², and filters by significance.  
- Top‑N export: writes your 20 best proxy ratios (by "r" and p < 0.05) to a CSV.

--------------------------------------------
USER OPTIONS
--------------------------------------------

These options are here for flexibility and ease of use and the purpose behind each is described in the code itself. In a nutshell, these allow you to control both the variables that the code will use to create the ratio, and the parameter that they are tested against. The parallel function is not necessary, but it does allow the code to run a lot faster. If you intend to test a large number of variables, it is advisable to have a reasonable max_k value. This limits the number of variables that can be used in combinations. If you do not set a max k value, then the number of combinations can quite easily reach into the trillions and would not be feasible to calculate unless you are using a supercomputer (which the parallel function would also run well on). To visually assist:


data  <-  state the name of your data file as it appears  # your abundance or proxy data
use_parallel   <- TRUE     # TRUE = parallel (no progress bar), FALSE = sequential (progress bar). This will make it faster on multi-core machines.
use_log_ratio  <- FALSE    # TRUE = log(A) - log(B), FALSE = A / B
max_k          <- 5        # maximum number of variables per numerator or denominator. Without a cap this can grow exponentially.
variable_cols  <- data[, x:y]   # columns (from x to y) in your data to consider for proxy combinations
environmental_col <- data$parameter_name  # the environmental parameter to correlate against, you can also do log(data$paramter_name) if that is convinient for you.
env_filter_min <- NA       # optional lower bound (set numeric to enable)
env_filter_max <- NA       # optional upper bound (set numeric to enable)
r_threshold    <- 0.60     # minimum |r| to include in results


Developed by Andrew Trott at ETH Zurich:
This project is licensed under the MIT License. See the LICENSE file for details.

