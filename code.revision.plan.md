Replace functions working with *stdDevSynthesisPrior to work with any of the prior related parameters.
- Examine code where string is used
- Possibly modify parameter object to include a prior category


## Usage in ROCParameter.cpp

```{cpp}
ROCParameter::ROCParameter(std::vector<double> stdDevSynthesisPrior, std::vector<unsigned> geneAssignment,
                                                std::vector<unsigned> _matrix, bool splitSer) : Parameter(22)

{
        unsigned _numMixtures = _matrix.size() / 2;
        std::vector<std::vector<unsigned>> thetaKMatrix;
        thetaKMatrix.resize(_numMixtures, std::vector<unsigned> (2, 0));
        unsigned index = 0;
        for (unsigned j = 0; j < 2; j++)
        {
                for (unsigned i = 0; i < _numMixtures; i++,index++)
                {
                        thetaKMatrix[i][j] = _matrix[index];
                }
        }
        initParameterSet(stdDevSynthesisPrior, _numMixtures, geneAssignment, thetaKMatrix, splitSer, "");
        initROCParameterSet();

}
```

It seems like we want to replace std::vector<double> stdDevSynthesisPrior with an object that holds the prior related parameters.
