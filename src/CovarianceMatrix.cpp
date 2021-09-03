#include "include/CovarianceMatrix.h"
#include "include/SequenceSummary.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif



//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


CovarianceMatrix::CovarianceMatrix()
{
    /* Initialize with numVariates = 2.
    // Equivalent to calling initCovarianceMatrix(2) */
    numVariates = 2;

    initCovarianceMatrix(numVariates);
}


CovarianceMatrix::CovarianceMatrix(unsigned _numVariates)
{
	numVariates = _numVariates;
	initCovarianceMatrix(_numVariates);
}


CovarianceMatrix::CovarianceMatrix(std::vector <double> &matrix)
{
    numVariates = (int)std::sqrt(matrix.size());
    covMatrix = matrix;
    choleskyMatrix.resize(matrix.size(), 0.0);
}


CovarianceMatrix::CovarianceMatrix(const CovarianceMatrix& other)
{
    numVariates = other.numVariates;
    covMatrix = other.covMatrix;
    choleskyMatrix = other.choleskyMatrix;
}


CovarianceMatrix& CovarianceMatrix::operator=(const CovarianceMatrix& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    numVariates = rhs.numVariates;
    covMatrix = rhs.covMatrix;
	choleskyMatrix = rhs.choleskyMatrix;
    return *this;
}


CovarianceMatrix& CovarianceMatrix::operator+(const CovarianceMatrix& rhs)
{
	std::vector<double> cov = rhs.covMatrix;
	for (unsigned i = 0; i < covMatrix.size(); i++)
	{
		covMatrix[i] += cov[i];
	}
	return *this;
}


CovarianceMatrix& CovarianceMatrix::operator*(const double &value)
{
	for (unsigned i = 0; i < covMatrix.size(); i++)
	{
		covMatrix[i] *= value;
	}
	return *this;
}


void CovarianceMatrix::operator*=(const double &value)
{
    for (unsigned i = 0; i < covMatrix.size(); i++)
    {
        covMatrix[i] *= value;
    }
}


bool CovarianceMatrix::operator==(const CovarianceMatrix& other) const 
{
    bool match = true;

    if (this->covMatrix != other.covMatrix) { match = false; }
    if (this->choleskyMatrix != other.choleskyMatrix) { match = false; }
    if (this->numVariates != other.numVariates) { match = false; }

    return match;
}


CovarianceMatrix::~CovarianceMatrix()
{
    //dtor
}





//--------------------------------------//
//---------- Matrix Functions ----------//
//--------------------------------------//


void CovarianceMatrix::initCovarianceMatrix(unsigned _numVariates)
{
    numVariates = _numVariates;
    unsigned vectorLength = numVariates * numVariates;
    covMatrix.resize(vectorLength);
    choleskyMatrix.resize(vectorLength);

	double diag_const = 0.01 / (double)numVariates;
    for (unsigned i = 0u; i < vectorLength; i++)
    {
        covMatrix[i] = (i % (numVariates + 1) ? 0.0 : diag_const);
        choleskyMatrix[i] = covMatrix[i];
    }
}


void CovarianceMatrix::setDiag(double val)
{
	for (unsigned i = 0u; i < covMatrix.size(); i++)
	{
		covMatrix[i] = (i % (numVariates + 1) ? covMatrix[i] : val);
	}
}


// adaptation of http://en.wikipedia.org/wiki/Cholesky_decomposition
// http://rosettacode.org/wiki/Cholesky_decomposition#C
void CovarianceMatrix::choleskyDecomposition()
{
    for (unsigned i = 0; i < numVariates; i++)
    {
        for (unsigned j = 0; j < (i + 1); j++)
        {
            double LsubstractSum = 0.0;
            for (unsigned k = 0; k < j; k++)
            {
                LsubstractSum += choleskyMatrix[i * numVariates + k] * choleskyMatrix[j * numVariates + k];
            }
            choleskyMatrix[i * numVariates + j] = (i == j) ? std::sqrt(covMatrix[i * numVariates + i] - LsubstractSum) :
                (1.0 / choleskyMatrix[j * numVariates + j]) * (covMatrix[i * numVariates + j] - LsubstractSum);
        }
    }
}


void CovarianceMatrix::printCovarianceMatrix()
{
    for (unsigned i = 0u; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0)
            my_print("\n");
        my_print("%\t", covMatrix[i]);
    }

    my_print("\n");
}


void CovarianceMatrix::printCholeskyMatrix()
{
    for (unsigned i = 0u; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0)
            my_print("\n");
        my_print("%\t", choleskyMatrix[i]);
    }

    my_print("\n");
}


std::vector<double>* CovarianceMatrix::getCovMatrix()
{
    std::vector<double> *ptr = &covMatrix;
    return ptr;
}


std::vector<double>* CovarianceMatrix::getCholeskyMatrix()
{
    std::vector<double> *ptr = &choleskyMatrix;
    return ptr;
}


int CovarianceMatrix::getNumVariates()
{
    return numVariates;
}


// This function multiplies the cholesky decomposition of the proposal matrix by the vector of indepdendent iid random variables (iidNumbers ~ N(0, 1))
// Thus it generates the shift in parameters from the current to the proposed values.
// Could be named proposedShiftCurrentParameters or the like
std::vector<double> CovarianceMatrix::transformIidNumbersIntoCovaryingNumbers(std::vector <double> iidNumbers)
{
    std::vector<double> covaryingNumbers;
    for (unsigned i = 0u; i < numVariates; i++)
    {
        double sum = 0.0;
        for (unsigned k = 0u; k < numVariates; k++)
        {
			// testing if [i * numVariates + k] or [k * numVariates + i], first option was default
            sum += choleskyMatrix[k * numVariates + i] * iidNumbers[k];
        }

        covaryingNumbers.push_back(sum);
    }
    return covaryingNumbers;
}


void CovarianceMatrix::calculateSampleCovariance(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned latestSample)
{
	//order of codonSpecificParameterTrace: paramType, category, numParam, samples
	unsigned numParamTypesInModel = (unsigned)codonSpecificParameterTrace.size();
	std::vector<unsigned> numCategoriesInModelPerParamType(numParamTypesInModel);
	// number of categories can vary between parameter types, see selection shared, mutation shared
	for (unsigned paramType = 0; paramType < numParamTypesInModel; paramType++)
	{
		numCategoriesInModelPerParamType[paramType] = (unsigned)codonSpecificParameterTrace[paramType].size();
	}


	unsigned start = latestSample - samples;
	
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);

	unsigned IDX = 0;
	for (unsigned paramType1 = 0; paramType1 < numParamTypesInModel; paramType1++)
	{
		unsigned numCategoriesInModel1 = numCategoriesInModelPerParamType[paramType1];
		for (unsigned category1 = 0; category1 < numCategoriesInModel1; category1++)
		{
			for (unsigned param1 = aaStart; param1 < aaEnd; param1++)
			{
				double mean1 = sampleMean(codonSpecificParameterTrace[paramType1][category1][param1], samples, latestSample);
				for (unsigned paramType2 = 0; paramType2 < numParamTypesInModel; paramType2++)
				{
					unsigned numCategoriesInModel2 = numCategoriesInModelPerParamType[paramType2];
					for (unsigned category2 = 0; category2 < numCategoriesInModel2; category2++)
					{
						for (unsigned param2 = aaStart; param2 < aaEnd; param2++)
						{
							double mean2 = sampleMean(codonSpecificParameterTrace[paramType2][category2][param2], samples, latestSample);
							double unscaledSampleCov = 0.0;
							for (unsigned i = start; i < latestSample; i++)
							{
								unscaledSampleCov += (codonSpecificParameterTrace[paramType1][category1][param1][i] - mean1) * (codonSpecificParameterTrace[paramType2][category2][param2][i] - mean2);
							}
							covMatrix[IDX] = unscaledSampleCov / ((double)samples - 1.0);

							IDX++;
						}
					}
				}
			}
		}
	}
}

void CovarianceMatrix::calculateSampleCovarianceForPANSE(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string codon, unsigned samples, unsigned latestSample)
{
    //order of codonSpecificParameterTrace: paramType, category, numParam, samples
    unsigned numParamTypesInModel = (unsigned)codonSpecificParameterTrace.size();
    std::vector<unsigned> numCategoriesInModelPerParamType(numParamTypesInModel);
    // number of categories can vary between parameter types, see selection shared, mutation shared
    for (unsigned paramType = 0; paramType < numParamTypesInModel; paramType++)
    {
        numCategoriesInModelPerParamType[paramType] = (unsigned)codonSpecificParameterTrace[paramType].size();
    }


    unsigned start = latestSample - samples;
    
    //unsigned aaStart, aaEnd;
    //SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
    unsigned codonIndex = SequenceSummary::codonToIndex(codon);
    unsigned IDX = 0;
    for (unsigned paramType1 = 0; paramType1 < numParamTypesInModel; paramType1++)
    {
        unsigned numCategoriesInModel1 = numCategoriesInModelPerParamType[paramType1];
        for (unsigned category1 = 0; category1 < numCategoriesInModel1; category1++)
        {
            double mean1 = sampleMean(codonSpecificParameterTrace[paramType1][category1][codonIndex], samples, latestSample,true);

            for (unsigned paramType2 = 0; paramType2 < numParamTypesInModel; paramType2++)
            {
                unsigned numCategoriesInModel2 = numCategoriesInModelPerParamType[paramType2];
                for (unsigned category2 = 0; category2 < numCategoriesInModel2; category2++)
                {               
                    double mean2 = sampleMean(codonSpecificParameterTrace[paramType2][category2][codonIndex], samples, latestSample,true);
                    double unscaledSampleCov = 0.0;
                    for (unsigned i = start; i < latestSample; i++)
                    {
                        unscaledSampleCov += (std::log(codonSpecificParameterTrace[paramType1][category1][codonIndex][i]) - mean1) * (std::log(codonSpecificParameterTrace[paramType2][category2][codonIndex][i]) - mean2);
                    }
                    covMatrix[IDX] = unscaledSampleCov / ((double)samples - 1.0);
                    IDX++;   
                }
            }
        }
    }
}


double CovarianceMatrix::sampleMean(std::vector<float> sampleVector, unsigned samples, unsigned latestSample,bool log_scale)
{
	double posteriorMean = 0.0;
	unsigned start = latestSample - samples;
	for (unsigned i = start; i < latestSample; i++)
	{
        if (log_scale)
        {
            posteriorMean += std::log(sampleVector[i]);
        }
        else
        {
		  posteriorMean += sampleVector[i];
	    }
    }
	return posteriorMean / (double)samples;
}

// Based on ramcmc package which is based on Vihola2012's method
// - Uses uniform iid variates u
// - Should be applied *after* the acceptance/rejection step
//
// ramcmc code
//    inline void adapt_S(arma::mat& S, arma::vec& u, double current, double target,
//      unsigned int n, double gamma) {
//    
//      double change = current - target;
//      u = S * u / arma::norm(u) * sqrt(std::min(1.0, u.n_elem * pow(n, -gamma)) *
//        std::abs(change));
//    
//      if(change > 0.0) {
//        chol_update(S, u);
//      } else {
//        chol_downdate(S, u);
//      }
//    }
//
//  
void CovarianceMatrix::adaptCholeskyMatrix(std::vector<double> &u, double current, double target, int updateNum)
{
  double error  = current - target;
  //Treat vector form of matrix as an array
  // See: https://stackoverflow.com/q/2923272/5322644
  std::vector<double> S = choleskyMatrix; // this compiles
  //std::vector<double> *S2 = &choleskyMatrix; // this compiles
  vector<double>::iterator i; 

  double uNorm = 0;

  for (i = u.begin(); i < u.end(); ++i)
  {
	  uNorm += u[i] * u[i];
  }
  
  uNorm = sqrt(uNorm);
  // ramcmc does the following
  // Why the u.n_elem? 
    //u = S * u/uNorm * sqrt(std::min(1.0, u.n_elem * pow(updateNum, -gamma)) * std::abs(error) );
  //first step of transformation
  double secondTerm = sqrt(std::min(1.0, u.n_elem * pow(updateNum, -gamma)) * std::abs(error) );
  for (i = u.begin(); i < u.end(); ++i)
  {
	  u[i] = u[i]/uNorm * secondTerm;
  }

 
  


  if(error > 0.0)  {
    updateCholesky(u);
  } else {
    downdateCholesky(u);
  }
}

// Based on ramcmc package which is based on Vihola2012's method
void checkCholeskyArgs()
{
}


void initCholesky() // may not need
{
  
}


// updateCholesky()
// - Based on ramcmc package which is based on Vihola2012's method
// - ramcmc code
//
//     inline arma::mat chol_update(arma::mat& L, arma::vec& u) {
//       unsigned int n = u.n_elem - 1;
//       for (arma::uword i = 0; i < n; i++) {             //uword = unsigned int
//         double r = sqrt(L(i,i) * L(i,i) + u(i) * u(i));
//         double c = r / L(i, i);
//         double s = u(i) / L(i, i);
//         L(i, i) = r;
//         L(arma::span(i + 1, n), i) =
//           (L(arma::span(i + 1, n), i) + s * u.rows(i + 1, n)) / c;
//         u.rows(i + 1, n) = c * u.rows(i + 1, n) -
//           s * L(arma::span(i + 1, n), i);
//       }
//       L(n, n) = sqrt(L(n, n) * L(n, n) + u(n) * u(n));
//       return L;
//     }
// Given the lower triangular matrix L obtained from the Cholesky decomposition of a 
// theoretical proposal matrix
// - ramcmc terminology
// - A is the proposal matrix
// - S is its Cholesky decomposition of A
// - L is a generic lower diagonal matrix
//   - That is, S is a specific instance of L
// - u is the vector of iid random variables, we use N(0,1)
// Updates choleskyMatrix (ramcmc uses L for general matrix to updateThis usually corresponds to the adapative S for proposal matrix) such that it corresponds to the decomposition of A + u*u'.
// Where u are the random iid numbers (iidNumbers in our code)  used to generate the proposed parameters (i.e. theta') or 
// Our code uses a flattened matrix (i.e. a 1 D vector) instead of a matrix 
void updateCholesky(std::vector<double> u) 
{
  //Treat vector form of matrix as an array
  // See: https://stackoverflow.com/q/2923272/5322644
  double *L = &choleskyMatrix.data(); 
  unsigned int n = numVariates - 1;//ramcmc: u.n_elem - 1;

  for (unsigned int i = 0u; i < n; i++) {
    double r = sqrt(L(i,i) * L(i,i) + u(i) * u(i)); // downdateCholesky uses "... - u(i) * u(i) 
    double c = r / L(i, i);
    double s = u(i) / L(i, i);
    L(i, i) = r;
    //    L(arma::span(i + 1, n), i) = (L(arma::span(i + 1, n), i) + s * u.rows(i + 1, n)) / c;
    for (unsigned j = i + 1; j <= n; j++) {
      L(j, i) = (L(j, i) + s * u(j)) / c;
      u(j) = c * u(j) - s * L(j, i);
    }
  }
  L(n, n) = sqrt(L(n, n) * L(n, n) + u(n) * u(n));
  return L;
}

//see updateCholesky for details
//
// ramcmc code
//   inline arma::mat chol_downdate(arma::mat& L, arma::vec& u) {
//      unsigned int n = u.n_elem - 1;
//      for (arma::uword i = 0; i < n; i++) {
//        double r = sqrt(L(i,i) * L(i,i) - u(i) * u(i));
//        double c = r / L(i, i);
//        double s = u(i) / L(i, i);
//        L(i, i) = r;
//        L(arma::span(i + 1, n), i) =
//          (L(arma::span(i + 1, n), i) - s * u.rows(i + 1, n)) / c;
//        u.rows(i + 1, n) = c * u.rows(i + 1, n) -
//          s * L(arma::span(i + 1, n), i);
//      }
//      L(n, n) = sqrt(L(n, n) * L(n, n) - u(n) * u(n));
//      return L;

void downdateCholesky(std::vector<double> u)
{
  //Treat vector form of matrix as an array
  // See: https://stackoverflow.com/q/2923272/5322644
  double *L = &choleskyMatrix.data(); 
  unsigned int n = numVariates - 1;//ramcmc: u.n_elem - 1;

  for (unsigned int i = 0u; i < n; i++) {
    double r = sqrt(L(i,i) * L(i,i) - u(i) * u(i)); // updateCholesky uses "... + u(i) * u(i) 
    double c = r / L(i, i);
    double s = u(i) / L(i, i);
    L(i, i) = r;
    //    L(arma::span(i + 1, n), i) = (L(arma::span(i + 1, n), i) + s * u.rows(i + 1, n)) / c;
    for (unsigned j = i + 1; j <= n; j++) {
      L(j, i) = (L(j, i) + s * u(j)) / c;
      u(j) = c * u(j) - s * L(j, i);
    }
  }
  L(n, n) = sqrt(L(n, n) * L(n, n) + u(n) * u(n));
  return L;
}


// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE

void CovarianceMatrix::setCovarianceMatrix(SEXP _matrix)
{
  std::vector<double> tmp;
  NumericMatrix matrix(_matrix);
  unsigned numRows = matrix.nrow();
  covMatrix.resize(numRows * numRows, 0.0);
  numVariates = numRows;
 
  //NumericMatrix stores the matrix by column, not by row. The loop
  //below transposes the matrix when it stores it.
  unsigned index = 0;
  for (unsigned i = 0; i < numRows; i++)
  {
    for (unsigned j = i; j < numRows * numRows; j += numRows, index++)
    {
      covMatrix[index] = matrix[j];
    }
  }
}





//----------------------------------//
//---------- RCPP Module -----------//
//----------------------------------//


RCPP_MODULE(CovarianceMatrix_mod)
{
  class_<CovarianceMatrix>( "CovarianceMatrix" )

        //Constructors & Destructors:
		.constructor("Empty Constructor")



		//Matrix Functions:
		.method("choleskyDecomposition", &CovarianceMatrix::choleskyDecomposition)
		.method("printCovarianceMatrix", &CovarianceMatrix::printCovarianceMatrix)
		.method("printCholeskyMatrix", &CovarianceMatrix::printCholeskyMatrix)
		.method("setCovarianceMatrix", &CovarianceMatrix::setCovarianceMatrix)
		;
}
#endif
