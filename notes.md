# Temporary file

- should move notes to labbook when appropriate.

Functions ROCModel::adapt*ProposalWidth in ROCModel.cpp seem to point to functions in Parameter.cpp

Implement Vihola's approach for 1D for sphi to start with 

- C++ Code
- adaptiveWidth vs. adaptationWidth: Synonyms or different?
  - Given functions that use adaptationWidth as a variable have adaptiveWidth in their name
- Is adaptiveWidth always scaled by thinning when it's set?
  -Right now this adjustment is done within the MCMCAlgorithm constructor
- variable lastIteration appears to actually be lastSample (or better, latestSample)
- Cusage: steps vs. iterations? Can we use one consistently? vs samples?
  - getSteps to adapt seem 
  - 

Selectively replaced via interactive vim sessions  with
    `for fname in $(git ls-files); do     grep -q "$PATTERN" "$fname" &&     vim -c "%s/$PATTERN/$REPLACEMENT/gc" -c 'wq' "$fname"; done`
    with 
    - PATTERN="sample"; REPLACEMENT="iteration"
    - PATTERN="sample steps"; REPLACEMENT="iteration"
    - PATTERN="sample steps"; REPLACEMENT="samples"
    - PATTERN="iteration (step)"; REPLACEMENT="iteration"
    - PATTERN="each step"; REPLACEMENT="each sample"
    - PATTERN="sample iterations to adapt"; REPLACEMENT="iterations to adapt"


What to do with getAdaptiveWidth which is in units of samples, not iterations


 /* getStepsToAdapt (RCPP EXPOSED)
 * Arguments: None
 * DEPRICATED FUNCTION REPLACED BY getIterationsToAdapt and getSamplesToAdapt
 * Return the value of iterationsToAdapt (formerly stepsToAdapt)
*/
//' @name getStepsToAdapt
//' @title getStepsToAdapt
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Return number of iterations (total iterations = samples * thinning) to allow proposal widths to adapt
//' @return number of sample steps to adapt
int MCMCAlgorithm::getStepsToAdapt()
{
	return iterationsToAdapt;
}


/* setStepsToAdapt (RCPP EXPOSED)
 * Arguments: iterations (unsigned)
 * DEPRICATED 
 * REPLACED BY setIterationsToAdapt and setSamplesToAdapt
 * Will set the specified iterations to adapt for the run if the value is less than samples * thinning (aka, the number
 * of iterations the run will last).The default parameter passed in as -1 uses the full iterations.
*/
//' @name setStepsToAdapt
//' @title setStepsToAdapt
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Set number of iterations (total iterations = samples * thinning) to allow proposal widths to adapt
//' @param steps a postive value
void MCMCAlgorithm::setStepsToAdapt(unsigned steps)
{
	if (steps <= samples * thinning)
		stepsToAdapt = steps;
	else
		my_printError("ERROR: Cannot set steps - value must be smaller than samples times thinning (maxIterations)\n");
}


/* getStepsToAdapt (RCPP EXPOSED)
 * Arguments: None
 * Return the value of stepsToAdapt
*/
//' @name getStepsToAdapt
//' @title getStepsToAdapt
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Return number of iterations (total iterations = samples * thinning) to allow proposal widths to adapt
//' @return number of sample steps to adapt
int MCMCAlgorithm::getStepsToAdapt()
{
	return iterationsToAdapt;
}



for fname in $(git ls-files); do
    grep -q 'PATTERN' "$fname" &&
    vim -c '%s/PATTERN/REPLACEMENT/gc' -c 'wq' "$fname"
done


- Replace setLastIteration with setLastSample and similar

- adaptiveWidth vs adaptationWidth


# 22 Jul 2021

- Why aren't the short Cov functions defined as inline functions?
  - Doing so means including 'inline' when explicitly defining the function with its code.
- In ramcmc/R/update.R,  the example uses `u <- runif(1, -1, 1)` which I find confusing.
  - That's because the are using unif(-s, s) as their proposal
    In ROCParameter.cpp we use 	`iidProposed.push_back(randNorm(0.0, 1.0))`
    - push_back() is a standard C++ function for vectors/arrays
    - C++ uses a `.` while R uses `$` for object oriented 
- ONE PULL REQUEST PER ISSUE!
- Updated testing
- src/testing.cpp is used to test C++ code
  - Alex believes this is called by testthat in R
- Issues
  - I should resolve
    - #364: "Move diag calculations"
    - #362: "Bug in AAToCodon"
      - Introduce exclude.reference arg
      - deprecate focal arg
  - Include issue number in branch name
    
- Currently p7utting update functions in CovarianceMatrix part of code.
  However, we want to be able to use the functions with individual, scaler parameters such as phi.
  Options
      - move functions to Parameter where they are applied to the covariance matrix and other parameters
      - Create an additional set of functions for the univariate case

# 26 Jul 2021

- CovarianceMatrix includes all categories
- numVariates in the sum of variates across all of the categories
   =  (numCodons * (numMutationCategories + numSelectionCategories)
- Alex's suggestion
  - Create vector in the parameter class to keep track of acceptance probabilities for each set of parameters
    - If the adaptiveWidth > step, use the average of these values instead of the single values as in most alogrithms
  - Within Cov object, have components for iterative updating of the cov matrix at each iteration
  - Update lambda only when proposal is being updated.
  - Proposal matrix is misnamed as `covariance matrix'
- Alex recommends stepping through code as it runs from R
- Make sure I use alpha = min(1, acceptance probability) when evaluating alpha
- For variable scope, check header files first
  - Parameter class
- sphi and fix_sphi is in Parameter.h
- fixDM is in ROCParameter.h and FONSEParameter.h
