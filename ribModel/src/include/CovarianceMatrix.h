#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <vector>
#ifndef STANDALONE
#include <Rcpp.h>
#endif
class CovarianceMatrix
{
    private:
    std::vector<double> covMatrix;
    std::vector<double> choleskiMatrix;
    int numVariates; //make static const again


    public:
        CovarianceMatrix();
				CovarianceMatrix(int _numVariates);
				CovarianceMatrix(std::vector <double> &matrix);
        virtual ~CovarianceMatrix();
				CovarianceMatrix(const CovarianceMatrix& other);
        CovarianceMatrix& operator=(const CovarianceMatrix& other);
        void initCovarianceMatrix(unsigned _numVariates);
        void choleskiDecomposition();
        void calculateCovarianceMatrixFromTraces(std::vector<std::vector <std::vector<double>>> trace, unsigned geneIndex, unsigned curSample, unsigned adaptiveWidte);
				void printCovarianceMatrix();
        void printCholeskiMatrix();
        void transformIidNumersIntoCovaryingNumbers(double* iidnumbers, double* covnumbers);
				#ifndef STANDALONE
				void setCovarianceMatrix(SEXP _matrix);
				#endif
    protected:
};

#endif // COVARIANCEMATRIX_H