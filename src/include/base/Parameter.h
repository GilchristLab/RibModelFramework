#ifndef PARAMETER_H
#define PARAMETER_H


#include "../Genome.h"
#include "../CovarianceMatrix.h"
#include "Trace.h"


#include <vector>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <ctime>
#include <sstream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

/* Note 1) -- on getSelectionCategory and getSynthesisRateCategory
 * These two functions are technically the same for readability.
 * Selection and Synthesis Rate are directly related even if they are not known
 * and thus are represented by the same variable. By splitting this
 * into Selection and Synthesis, avoids confusing the two, however.
*/

class Parameter {
	private:

		//STATICS - Sorting Functions:
		void quickSortPair(double a[], int b[], int first, int last);
		static int pivotPair(double a[], int b[], int first, int last);

		

		std::vector<double> codonSpecificPrior;

		bool fix_stdDevSynthesis = false;
	
	public:

		static const std::string allUnique;
		static const std::string selectionShared;
		static const std::string mutationShared;
		static const std::string elongationShared;
		static const std::string nseShared;

		static const unsigned dM;
		static const unsigned dEta;
		static const unsigned dOmega;
		static const unsigned alp;
		static const unsigned lmPri;
		static const unsigned nse;

		// Phi prior selection codes (task #12: mixture phi prior).
		// Defaults map to existing behavior (SINGLE_LN + MEAN).
		static const unsigned PHI_PRIOR_SINGLE_LN;   // = 0; legacy single LogNormal (default)
		static const unsigned PHI_PRIOR_MIXTURE_LN;  // = 1; mixture of two LogNormals
		static const unsigned PHI_CONSTRAINT_MEAN;    // = 0; anchors E[phi]   = 1 (default)
		static const unsigned PHI_CONSTRAINT_MEDIAN;  // = 1; anchors median[phi] = 1

#ifdef STANDALONE
		static std::default_random_engine generator; // static to make sure that the same generator is used during the runtime.
#endif


		//Constructors & Destructors:
		Parameter();
		Parameter(unsigned maxGrouping);
		Parameter& operator=(const Parameter& rhs);
		virtual ~Parameter();


		//Initialization and Restart Functions: TODO: test
		void initParameterSet(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures,
			std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> mixtureDefinitionMatrix,
			bool splitSer = true, std::string _mutationSelectionState = "allUnique"); //Mostly tested; TODO caveats
		void initBaseValuesFromFile(std::string filename);
		void writeBasicRestartFile(std::string filename);
		void initCategoryDefinitions(std::string mutationSelectionState,
			std::vector<std::vector<unsigned>> mixtureDefinitionMatrix);
		void InitializeSynthesisRate(Genome& genome, double sd_phi);
		void InitializeSynthesisRate(double sd_phi);
		void InitializeSynthesisRate(std::vector<double> expression);
		std::vector<double> readPhiValues(std::string filename); //General function, possibly move


		//Prior functions: TODO: test
		double getCodonSpecificPriorStdDev(unsigned paramType);


		//Mixture Definition Matrix and Category Functions: Mostly tested, see comments
		void setNumMutationSelectionValues(std::string mutationSelectionState,
			std::vector<std::vector<unsigned>> mixtureDefinitionMatrix); //TODO: test
		void printMixtureDefinitionMatrix(); //Untested
		double getCategoryProbability(unsigned mixtureElement);
		void setCategoryProbability(unsigned mixtureElement, double value);
		unsigned getNumMutationCategories(); //TODO caveat
		unsigned getNumSelectionCategories(); //TODO caveat
		unsigned getNumSynthesisRateCategories(); //TODO caveat
		unsigned getMutationCategory(unsigned mixtureElement);
		unsigned getSelectionCategory(unsigned mixtureElement); //see Note 1) at top of file.
		unsigned getSynthesisRateCategory(unsigned mixtureElement); //see Note 1) at top of file.
		std::vector<unsigned> getMixtureElementsOfMutationCategory(unsigned category); //TODO caveat
		std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned category); //TODO caveat
		std::vector<unsigned> getMixtureElementsOfSynthesisRateCategory(unsigned category); //TODO caveat
		std::string getMutationSelectionState();
		unsigned getNumAcceptForCspForIndex(unsigned i); //Only for unit testing.


		//Group List Functions: All tested
		void setGroupList(std::vector<std::string> gl);
		std::string getGrouping(unsigned index);
		std::vector<std::string> getGroupList();
		unsigned getGroupListSize();


		//stdDevSynthesisRate Functions: Mostly tested, see comments.
		double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false);
		virtual void proposeStdDevSynthesisRate(); //TODO: test
		void fixStdDevSynthesis();
		void setStdDevSynthesisRate(double stdDevSynthesisRate, unsigned selectionCategory);
		double getCurrentStdDevSynthesisRateProposalWidth();
		unsigned getNumAcceptForStdDevSynthesisRate(); //Only for unit testing.
		void updateStdDevSynthesisRate(); //TODO: test
		double getStdCspForIndex(unsigned i); //Only for unit testing.


		// Phi prior type / mixture-LN storage (task #12a). Defaults preserve
		// legacy behavior (single LogNormal anchored at E[phi]=1).
		unsigned getPhiPriorType();
		void setPhiPriorType(unsigned type);
		unsigned getPhiPriorConstraint();
		void setPhiPriorConstraint(unsigned constraint);

		// Per-mixture-category mixture-LN parameters.
		double getPhiMixtureP(unsigned mixtureCategory, bool proposed = false);
		void setPhiMixtureP(double p, unsigned mixtureCategory);
		double getPhiMixtureMu1(unsigned mixtureCategory, bool proposed = false);
		void setPhiMixtureMu1(double mu1, unsigned mixtureCategory);
		double getPhiMixtureSigma1(unsigned mixtureCategory, bool proposed = false);
		void setPhiMixtureSigma1(double sigma1, unsigned mixtureCategory);
		double getPhiMixtureSigma2(unsigned mixtureCategory, bool proposed = false);
		void setPhiMixtureSigma2(double sigma2, unsigned mixtureCategory);
		// mu2 is not stored; derived from the constraint at evaluation time.
		double getPhiMixtureMu2Derived(unsigned mixtureCategory);

		// Hyperparameter accessors.
		double getPhiMixtureHyperPAlpha();
		void setPhiMixtureHyperPAlpha(double v);
		double getPhiMixtureHyperPBeta();
		void setPhiMixtureHyperPBeta(double v);
		double getPhiMixtureHyperMu1Mean();
		void setPhiMixtureHyperMu1Mean(double v);
		double getPhiMixtureHyperMu1Sd();
		void setPhiMixtureHyperMu1Sd(double v);
		double getPhiMixtureHyperSigma1Scale();
		void setPhiMixtureHyperSigma1Scale(double v);
		double getPhiMixtureHyperSigma2Scale();
		void setPhiMixtureHyperSigma2Scale(double v);

		// Resize storage to numSynthesisRateCategories and seed defaults.
		// Called from initParameterSet; safe to call again to re-seed.
		void initPhiMixtureStorage();

		// Log-prior on phi at the given mixture category, switching internally
		// on phiPriorType. Single source of truth for the prior choice -- both
		// MCMCAlgorithm and per-model codon updates route through this.
		// For SINGLE_LN this is bit-for-bit identical to the legacy inline
		// LogNormal(-sigma^2/2, sigma) computation.
		double getLogPhiPrior(double phi, unsigned mixtureCategory);

		// Single-site Metropolis-Hastings update for the mixture hyperparams
		// (p, mu1, sigma1, sigma2), per mixture category (task #12c.1).
		// Sequential transformed-scale random walks with Jacobian corrections:
		//   p     -> logit-scale (jacobian: log(p*(1-p)))
		//   mu1   -> identity-scale
		//   sigma -> log-scale (jacobian: log(sigma))
		// No-op when phiPriorType != PHI_PRIOR_MIXTURE_LN. Fixed proposal
		// widths (std_phiMixture*); adaptive tuning + trace storage land
		// in 12c.2. Mirrors prototypes/phi_mixture.R::mh_phi_mixture.
		void updatePhiMixtureHyperparameters(Genome& genome);

		// Accept counters for the four mixture hyperparams. Reset each
		// adaptive-width window (task #12c.2).
		unsigned getNumAcceptForPhiMixtureP();
		unsigned getNumAcceptForPhiMixtureMu1();
		unsigned getNumAcceptForPhiMixtureSigma1();
		unsigned getNumAcceptForPhiMixtureSigma2();

		// Write the current (p, mu1, sigma1, sigma2) per mixture category into
		// the trace at this sample index. No-op for SINGLE_LN (task #12c.2).
		void updatePhiMixtureTrace(unsigned sample);

		// Adapt the four proposal widths toward target acceptance ~0.25, mirror
		// of adaptStdDevSynthesisRateProposalWidth. Pushes one acceptance-rate
		// trace point per param and resets accept counters (task #12c.2).
		void adaptPhiMixtureProposalWidths(unsigned adaptationWidth, bool adapt);

		// Returns the gene count implied by the parameter object's internal
		// state (the inner size of currentSynthesisRateLevel).  Used by
		// MCMCAlgorithm::run to validate that the genome it is about to
		// iterate over has the same number of genes the parameter object
		// was sized for; mismatch was previously a silent OOB segfault.
		// Returns 0 for parameter objects that have not been sized yet
		// (e.g. fresh-constructed by the default ctor before initParameterSet).
		unsigned getNumGenesFromState() const {
			return currentSynthesisRateLevel.empty()
				? 0u
				: (unsigned) currentSynthesisRateLevel[0].size();
		}

		// Restart-file build-info accessors.  Captures provenance of the
		// .rst that produced this parameter object (version, commit SHA,
		// configure-time build date, and the timestamp the .rst was
		// written).  Empty strings (and generation LEGACY_2022_RELEASE /
		// MODERN_NO_BUILDINFO / UNKNOWN) when the .rst lacked a
		// >buildInfo: block.  See docs/VERSIONING.md for the schema.
		enum RestartFileGeneration {
			REST_GEN_UNKNOWN              = 0,
			REST_GEN_LEGACY_2022_RELEASE  = 1,
			REST_GEN_MODERN_NO_BUILDINFO  = 2,
			REST_GEN_MODERN_WITH_BUILDINFO = 3
		};
		std::string getRestartFileVersion()    const { return restartFileVersion; }
		std::string getRestartFileCommitSha()  const { return restartFileCommitSha; }
		std::string getRestartFileBuildDate()  const { return restartFileBuildDate; }
		std::string getRestartFileWrittenAt()  const { return restartFileWrittenAt; }
		unsigned    getRestartFileGeneration() const { return (unsigned) restartFileGeneration; }
		std::string getRestartFileGenerationName() const;


		//Synthesis Rate Functions: Mostly tested, see comments
		double getSynthesisRate(unsigned geneIndex, unsigned mixtureElement, bool proposed = false);
		double getCurrentSynthesisRateProposalWidth(unsigned expressionCategory, unsigned geneIndex);
		double getSynthesisRateProposalWidth(unsigned geneIndex, unsigned mixtureElement);
		void proposeSynthesisRateLevels(); //TODO: test
		void setSynthesisRate(double phi, unsigned geneIndex, unsigned mixtureElement);
		void updateSynthesisRate(unsigned geneIndex); //TODO: test
		void updateSynthesisRate(unsigned geneIndex, unsigned mixtureElement); //TODO: test
		unsigned getNumAcceptForSynthesisRate(unsigned expressionCategory, unsigned geneIndex); //Only for unit testing


		//Noise Functions...updating AnaCoDa to allow all models (instead of just ROC) to use empirical gene expression values to inform estimation of \phi

		double getObservedSynthesisNoise(unsigned index);
		void setObservedSynthesisNoise(unsigned index, double se);

		//noiseOffset Functions:
		double getNoiseOffset(unsigned index, bool proposed = false);
		double getCurrentNoiseOffsetProposalWidth(unsigned index);
		void proposeNoiseOffset();
		void setNoiseOffset(unsigned index, double _NoiseOffset);
		void updateNoiseOffset(unsigned index);
		void updateGibbsSampledHyperParameters(Genome &genome, bool withPhi,bool fix_sEpsilon);

		void setNumObservedPhiSets(unsigned _phiGroupings);

		// noise Functions:
		void setInitialValuesForSepsilon(std::vector<double> seps);

		//Posterior, Variance, and Estimates Functions for noise:
		double getNoiseOffsetPosteriorMean(unsigned index, unsigned samples);
		double getNoiseOffsetVariance(unsigned index, unsigned samples, bool unbiased = true);

		//Adaptive Width Functions:
		void adaptNoiseOffsetProposalWidth(unsigned adaptationWidth, bool adapt);

		//Iteration Functions: All tested
		unsigned getLastIteration();
		void setLastIteration(unsigned iteration);


		//Trace Functions: TODO: test
		Trace& getTraceObject();
		void setTraceObject(Trace _trace);
		void updateObservedSynthesisNoiseTraces(unsigned sample);
		void updateNoiseOffsetTraces(unsigned sample);
		void updateStdDevSynthesisRateTrace(unsigned sample);
		void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex);
		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex);
		void updateMixtureProbabilitiesTrace(unsigned samples);


		//Adaptive Width Functions: TODO: test
		void adaptStdDevSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt);
		void adaptSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt);


		//Posterior, Variance, and Estimates Functions: TODO: test
		double getStdDevSynthesisRatePosteriorMean(unsigned samples, unsigned mixture);
		double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, bool log_scale=false);

		double getCodonSpecificPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon,
			unsigned paramType, bool withoutReference = true, bool byGene = false, bool log_scale = false);
		double getStdDevSynthesisRateVariance(unsigned samples, unsigned mixture, bool unbiased);
		double getSynthesisRateVariance(unsigned samples, unsigned geneIndex,
			bool unbiased = true, bool log_scale = false);
		double getCodonSpecificVariance(unsigned mixtureElement, unsigned samples, std::string &codon,
			unsigned paramType, bool unbiased, bool withoutReference = true, bool log_scale = false);
	        std::vector<double> getCodonSpecificQuantile(unsigned mixtureElement, unsigned samples, std::string &codon,
			unsigned paramType, std::vector<double> probs, bool withoutReference, bool log_scale = false);
		std::vector<double> getExpressionQuantile(unsigned samples, unsigned geneIndex,
			std::vector<double> probs, bool log_scale = false);
		std::vector<double> calculateQuantile(std::vector<float> &parameterTrace, unsigned samples, std::vector<double> probs, bool log_scale=false);
		unsigned getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex);
		std::vector<double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);


		//Other Functions: Mostly tested, see comments
		unsigned getNumParam();
		unsigned getNumMixtureElements();
		unsigned getNumElongationMixtures();
		unsigned getNumObservedPhiSets();
		void setMixtureAssignment(unsigned gene, unsigned value);
		unsigned getMixtureAssignment(unsigned gene);
		virtual std::vector <std::vector <double> > calculateSelectionCoefficients(unsigned sample); //TODO: test



		//Static Functions: TODO: test
		static double calculateSCUO(Gene& gene);
		static void drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b),
			double* randomNumbers);
		static void drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r), double* randomNumber);
		static double randNorm(double mean, double sd);
		static double randLogNorm(double m, double s);
		static double randExp(double r);
		static double randGamma(double shape, double rate);
		static void randDirichlet(std::vector <double> &input, unsigned numElements, std::vector <double> &output);
		static double randUnif(double minVal, double maxVal);
		static unsigned randMultinom(double *probabilities, unsigned mixtureElements);
		static unsigned randMultinom(std::vector <double> &probabilities, unsigned mixtureElements);
		static double densityNorm(double x, double mean, double sd, bool log = false);
		static double densityLogNorm(double x, double mean, double sd, bool log = false);

		// Closed-form derivation of mu2 for the mixture phi prior under either
		// the mean=1 or median=1 constraint. See prototypes/phi_mixture.R for
		// the math derivation. Returns NaN if the constraint is infeasible
		// (the lower component alone already overshoots the anchor target).
		// constraint: PHI_CONSTRAINT_MEAN or PHI_CONSTRAINT_MEDIAN.
		static double deriveMu2(double p, double mu1, double sigma1, double sigma2,
		                        unsigned constraint);

		// Log-density of the constrained mixture-LN at x. Returns -DBL_MAX
		// (the "impossible" sentinel used by densityLogNorm) when the
		// derived mu2 is infeasible or violates the label-switching guard
		// mu2 >= mu1. M-H ratios using this sentinel auto-reject the move.
		static double densityLogNormMixture(double x, double p, double mu1,
		                                    double sigma1, double sigma2,
		                                    unsigned constraint, bool log = false);
		//double getMixtureAssignmentPosteriorMean(unsigned samples, unsigned geneIndex);
		// TODO: implement variance function, fix Mean function (won't work with 3 groups)





		//R Section:

#ifndef STANDALONE

		//Initialization and Restart Functions:
		void initializeSynthesisRateByGenome(Genome& genome, double sd_phi);
		void initializeSynthesisRateByRandom(double sd_phi);
		void initializeSynthesisRateByList(std::vector<double> expression);
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);



		//Mixture Definition Matrix and Category Functions:
		unsigned getMutationCategoryForMixture(unsigned mixtureElement);
		unsigned getSelectionCategoryForMixture(unsigned mixtureElement);
		unsigned getSynthesisRateCategoryForMixture(unsigned mixtureElement);
		std::vector<std::vector<unsigned>> getCategories();
		void setCategories(std::vector<std::vector<unsigned>> _categories);
		void setCategoriesForTrace();
		void setNumMutationCategories(unsigned _numMutationCategories);
		void setNumSelectionCategories(unsigned _numSelectionCategories);



		//Synthesis Rate Functions:
		std::vector<std::vector<double>> getSynthesisRateR();
		std::vector<double> getCurrentSynthesisRateForMixture(unsigned mixture);



		//Posterior, Variance, and Estimates Functions:
		double getSynthesisRatePosteriorMeanForGene(unsigned samples, unsigned geneIndex, bool log_scale);
		double getSynthesisRateVarianceForGene(unsigned samples, unsigned geneIndex, bool unbiased, bool log_scale);
		unsigned getEstimatedMixtureAssignmentForGene(unsigned samples, unsigned geneIndex);

		std::vector<double> getEstimatedMixtureAssignmentProbabilitiesForGene(unsigned samples, unsigned geneIndex);

		double getCodonSpecificPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon,
			unsigned paramType, bool withoutReference, bool log_scale = false);
		double getCodonSpecificVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon,
			unsigned paramType, bool unbiased, bool withoutReference, bool log_scale = false);
        	std::vector<double> getCodonSpecificQuantileForCodon(unsigned mixtureElement, unsigned samples,
        		std::string &codon, unsigned paramType, std::vector<double> probs, bool withoutReference, bool log_scale = false);
		std::vector<double> getExpressionQuantileForGene(unsigned samples,
			unsigned geneIndex, std::vector<double> probs, bool log_scale);



		//Other Functions:
		SEXP calculateSelectionCoefficientsR(unsigned sample);
		std::vector<unsigned> getMixtureAssignmentR();
		void setMixtureAssignmentR(std::vector<unsigned> _mixtureAssignment);
		unsigned getMixtureAssignmentForGeneR(unsigned geneIndex);
		void setMixtureAssignmentForGene(unsigned geneIndex, unsigned value);
		void setNumMixtureElements(unsigned _numMixtures);
		void setNumElongationMixtures(unsigned _numElongationMixtures);

#endif

	protected:
		Trace traces;

		unsigned adaptiveStepPrev;
		unsigned adaptiveStepCurr;
		
		std::vector<CovarianceMatrix> covarianceMatrix;
		std::vector<mixtureDefinition> categories;
		std::vector<double> categoryProbabilities;
		std::vector<std::vector<unsigned>> mutationIsInMixture;
		std::vector<std::vector<unsigned>> selectionIsInMixture;
		std::vector<std::vector<unsigned>> phiIsInMixture;
		unsigned numMutationCategories; //TODO Probably needs to be renamed
		unsigned numSelectionCategories; //TODO Probably needs to be renamed
		unsigned numSynthesisRateCategories;


		std::vector<unsigned> numAcceptForCodonSpecificParameters;
		std::string mutationSelectionState; //TODO: Probably needs to be renamed

        //<Alpha or Lambda or Mutation or Selection < Mixture < Codon >>>
		std::vector<std::vector<std::vector<double>>> proposedCodonSpecificParameter;
		std::vector<std::vector<std::vector<double>>> currentCodonSpecificParameter;

		std::vector<unsigned> mixtureAssignment;
		std::vector<std::string> groupList;
		unsigned maxGrouping;


		std::vector<double> stdDevSynthesisRate_proposed;
		std::vector<double> stdDevSynthesisRate;
		double bias_stdDevSynthesisRate; //NOTE: Currently, this value is always set to 0.0
		double std_stdDevSynthesisRate;
		unsigned numAcceptForStdDevSynthesisRate;
		std::vector<double> std_csp;

		// Phi prior selection (task #12). Default values map to legacy behavior:
		// SINGLE_LN + MEAN constraint == LogNormal(-sigma^2/2, sigma).
		unsigned phiPriorType;
		unsigned phiPriorConstraint;

		// Mixture-LN phi prior, per-mixture-category storage. Empty unless
		// phiPriorType == PHI_PRIOR_MIXTURE_LN; resized in initPhiMixtureStorage().
		// mu2 is NOT stored: it is derived from (p, mu1, sigma1, sigma2, constraint)
		// at evaluation time.
		std::vector<double> phiMixtureP;
		std::vector<double> phiMixtureP_proposed;
		std::vector<double> phiMixtureMu1;
		std::vector<double> phiMixtureMu1_proposed;
		std::vector<double> phiMixtureSigma1;
		std::vector<double> phiMixtureSigma1_proposed;
		std::vector<double> phiMixtureSigma2;
		std::vector<double> phiMixtureSigma2_proposed;

		// Per-param proposal widths and accept counters (used in 12c).
		double std_phiMixtureP;
		double std_phiMixtureMu1;
		double std_phiMixtureSigma1;
		double std_phiMixtureSigma2;
		unsigned numAcceptForPhiMixtureP;
		unsigned numAcceptForPhiMixtureMu1;
		unsigned numAcceptForPhiMixtureSigma1;
		unsigned numAcceptForPhiMixtureSigma2;

		// Hyperparameters (scalar; settable from R). Defaults match the
		// validated R prototype (Beta(8,2), Normal(0,10), half-Normal(1)).
		double phiMixtureHyper_p_alpha;
		double phiMixtureHyper_p_beta;
		double phiMixtureHyper_mu1_mean;
		double phiMixtureHyper_mu1_sd;
		double phiMixtureHyper_sigma1_scale;
		double phiMixtureHyper_sigma2_scale;


		std::vector <double> observedSynthesisNoise;

		std::vector <double> noiseOffset_proposed;
		std::vector <double> noiseOffset; //A_Phi
		std::vector <double> std_NoiseOffset;
		std::vector <double> numAcceptForNoiseOffset;
		
        //Unknown indexing hoping (mixture) then gene
		std::vector<std::vector<double>> proposedSynthesisRateLevel;
		std::vector<std::vector<double>> currentSynthesisRateLevel;
		std::vector<std::vector<unsigned>> numAcceptForSynthesisRate;

		unsigned lastIteration;

		unsigned int numParam;
		unsigned numMixtures;
		unsigned numElongationMixtures; // Intended for just PA and PANSE for now, but can begin thinking of expanding to other models.
		unsigned obsPhiSets;

		double bias_phi; //NOTE: Currently, this value is always set to 0.0
		std::vector<std::vector<double>> std_phi;

		// Restart-file build-info read from the >buildInfo: block at the
		// top of the file.  Empty when the file lacked the block (e.g.
		// the official 2022 release, which writes neither >buildInfo:
		// nor >numSynthesisRateCategories:).  See Parameter::detectRestartFileGeneration().
		std::string restartFileVersion;
		std::string restartFileCommitSha;
		std::string restartFileBuildDate;
		std::string restartFileWrittenAt;
		RestartFileGeneration restartFileGeneration;

};

#endif // PARAMETER_H
