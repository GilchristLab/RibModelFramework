#ifndef ROCMODEL_H
#define ROCMODEL_H

#include "../base/Model.h"
#include "ROCParameter.h"

class ROCModel : public Model
{
    private:
		ROCParameter *parameter;

		double calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue);
		double calculateMutationPrior(std::string grouping, bool proposed = false); // TODO add to FONSE as well? // cedric
		double calculateSelectionPrior(std::string grouping, bool proposed = false); 
		void obtainCodonCount(SequenceSummary *sequenceSummary, std::string curAA, int codonCount[]);
    public:
		//Constructors & Destructors:
		ROCModel(bool _withPhi = false, bool _fix_sEpsilon = false);
		virtual ~ROCModel();

		std::string type = "ROC";

		//Likelihood Ratio Functions:
		virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k,
					double* logProbabilityRatio);
		virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
					std::vector<double> &logAcceptanceRatioForAllMixtures, std::string param="Evolutionary");
		virtual void calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration,
					std::vector <double> &logProbabilityRatio);

		// Full L(data | theta) evaluation at the model's CURRENT parameter
		// state.  Mirrors PA/PANSE Model::calculateLogLikelihood signature.
		// Sums per-AA per-gene contributions across the entire genome via
		// calculateLogLikelihoodPerAAPerGene; reads dM, dEta, and phi from
		// the bound parameter object rather than taking them as arguments.
		//
		// For DIC: call once at the posterior-mean parameter state to get
		//          D(theta_bar); D_bar comes for free from the existing
		//          MCMC$getLogLikelihoodTrace().
		// For bridge sampling: call once per proposal-distribution sample.
		//
		// Currently assumes numMixtures == 1 (the Lokiarchaeota case);
		// multi-mixture marginalization would require summing across
		// mixture categories weighted by mixture probabilities, which
		// the MCMC accept-reject pathway handles internally but is not
		// yet abstracted here.
		//
		// LIMITATION: Works on a Parameter that's been freshly initialized
		// via initializeParameterObject() (in-memory, post-MCMC).  Does
		// NOT work on a Parameter loaded from .Rdata via
		// loadParameterObject(): writeParameterObject() does not preserve
		// currentSynthesisRateLevel, so the per-gene phi accessor returns
		// uninitialized state.  For post-hoc analysis from saved .Rdata,
		// either (a) restore phi from trace via R helper, or (b) call this
		// method from within the original R session that produced the fit.
		double calculateLogLikelihood(Genome& genome);


		//Initialization and Restart Functions:
		virtual void initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate = true);
		virtual void writeRestartFile(std::string filename);



		//Category Functions:
		virtual double getCategoryProbability(unsigned i);
		virtual unsigned getMutationCategory(unsigned mixture);
		virtual unsigned getSelectionCategory(unsigned mixture);
		virtual unsigned getSynthesisRateCategory(unsigned mixture);
		virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k);



		//Group List Functions:
		virtual unsigned getGroupListSize();
		virtual std::string getGrouping(unsigned index);



		//stdDevSynthesisRate Functions:
		virtual double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false);
		virtual double getCurrentStdDevSynthesisRateProposalWidth();
		virtual void updateStdDevSynthesisRate();
		virtual double getLogPhiPrior(double phi, unsigned mixtureCategory);
		virtual void updatePhiMixtureHyperparameters(Genome& genome);



		//Synthesis Rate Functions:
		virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false);
		virtual void updateSynthesisRate(unsigned i, unsigned k);



		//Iteration Functions:
		virtual unsigned getLastIteration();
		virtual void setLastIteration(unsigned iteration);



		//Trace Functions:
		virtual void updateStdDevSynthesisRateTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i) ;
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i) ;
		virtual void updateMixtureProbabilitiesTrace(unsigned sample);
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping);
		virtual void updateHyperParameterTraces(unsigned sample);
		virtual void updateTracesWithInitialValues(Genome &genome);



		//Adaptive Width Functions:
		virtual void adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastSample, bool adapt = true);
		virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt = true);



		//Other Functions:
		virtual void proposeCodonSpecificParameter();
		virtual void proposeHyperParameters();
		virtual void proposeSynthesisRateLevels();

		virtual unsigned getNumPhiGroupings() ;
		virtual unsigned getMixtureAssignment(unsigned index);
		virtual unsigned getNumMixtureElements() ;
		virtual unsigned getNumSynthesisRateCategories();

		virtual void setNumPhiGroupings(unsigned value);
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene);
		virtual void setCategoryProbability(unsigned mixture, double value);

		virtual void updateCodonSpecificParameter(std::string grouping);
		virtual void updateCodonSpecificParameter(std::string grouping, std::string param = "Evolutionary");

		//virtual void updateGibbsSampledHyperParameters(Genome &genome);
		virtual void updateAllHyperParameter();
		virtual void updateHyperParameter(unsigned hp);

		void simulateGenome(Genome &genome);
		virtual void printHyperParameters();
		ROCParameter getParameter();
		void setParameter(ROCParameter &_parameter);

		// Override: delegate to the typed parameter pointer (which shadows
		// Model::parameter).  See base class docs for rationale.
		virtual unsigned getNumGenesFromState() const {
			return parameter ? parameter->getNumGenesFromState() : 0u;
		}
		virtual double calculateAllPriors(bool proposed=false);
		void calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[]);
		void calculateLogCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[]);
		virtual void getParameterForCategory(unsigned category, unsigned param, std::string aa, bool proposal, double* returnValue);


	
		 double getNoiseOffset(unsigned index, bool proposed = false);
		 double getObservedSynthesisNoise(unsigned index) ;
		 double getCurrentNoiseOffsetProposalWidth(unsigned index);
		 void updateNoiseOffset(unsigned index);
		 void updateNoiseOffsetTrace(unsigned sample);
		 void updateObservedSynthesisNoiseTrace(unsigned sample);
		 void adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		 void updateGibbsSampledHyperParameters(Genome &genome);

		virtual bool getParameterTypeFixed(std::string csp_parameters);
		virtual bool isShared(std::string csp_parameters);
		// Override base Model::recordCSPStepAlpha so it reaches ROCModel's
		// own typed parameter ptr (the base Model::parameter is shadowed
		// and stays nullptr).  Used by Vihola2012CSPAdapter for per-step
		// alpha capture.
		virtual void recordCSPStepAlpha(std::string grouping, double alpha) override;


		//R Section:
#ifndef STANDALONE
		std::vector<double> CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection,
							double phi);
#endif //STANDALONE

    protected:
    	
};

#endif // ROCMODEL_H
