
#ifndef STANDALONE
#include "include/base/Trace.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_MODULE(Trace_mod)
{
  class_<Trace>( "Trace" )
  
    //Getter Functions:
    .method("getStdDevSynthesisRateAcceptanceRateTrace", &Trace::getStdDevSynthesisRateAcceptanceRateTrace)
    .method("getSynthesisRateTrace", &Trace::getSynthesisRateTrace)
    .method("getSynthesisRateAcceptanceRateTrace", &Trace::getSynthesisRateAcceptanceRateTrace)
    .method("getCodonSpecificAcceptanceRateTraceForAA", &Trace::getCodonSpecificAcceptanceRateTraceForAA)
    .method("getCodonSpecificAcceptanceRateTraceForCodon", &Trace::getCodonSpecificAcceptanceRateTraceForCodon)
    .method("getMixtureAssignmentTrace", &Trace::getMixtureAssignmentTrace)
    .method("getCodonSpecificAcceptanceRateTrace", &Trace::getCodonSpecificAcceptanceRateTrace)
    .method("getNseRateSpecificAcceptanceRateTrace", &Trace::getNseRateSpecificAcceptanceRateTrace)
    .method("getMixtureProbabilitiesTrace", &Trace::getMixtureProbabilitiesTrace)
    .method("getExpectedSynthesisRateTrace", &Trace::getExpectedSynthesisRateTrace)
    .method("getSynthesisOffsetAcceptanceRateTrace", &Trace::getSynthesisOffsetAcceptanceRateTrace)
    .method("getSynthesisOffsetAcceptanceRateTraceForIndex", &Trace::getSynthesisOffsetAcceptanceRateTraceForIndex)
    .method("getCodonSpecificParameterTrace", &Trace::getCodonSpecificParameterTraceByParamType)
    .method("getSynthesisRateAcceptanceRateTraceByMixtureElementForGene",
            &Trace::getSynthesisRateAcceptanceRateTraceByMixtureElementForGeneR)
    .method("getSynthesisRateTraceForGene", &Trace::getSynthesisRateTraceForGeneR)
    .method("getSynthesisRateTraceByMixtureElementForGene", &Trace::getSynthesisRateTraceByMixtureElementForGeneR)
    .method("getMixtureAssignmentTraceForGene", &Trace::getMixtureAssignmentTraceForGeneR)
    .method("getMixtureProbabilitiesTraceForMixture", &Trace::getMixtureProbabilitiesTraceForMixtureR)
    .method("getStdDevSynthesisRateTraces", &Trace::getStdDevSynthesisRateTraces)
    .method("getNumberOfMixtures", &Trace::getNumberOfMixtures)
    // Phi mixture hyperparam traces (task #12c.2)
    .method("getPhiMixturePTrace",      &Trace::getPhiMixturePTrace)
    .method("getPhiMixtureMu1Trace",    &Trace::getPhiMixtureMu1Trace)
    .method("getPhiMixtureSigma1Trace", &Trace::getPhiMixtureSigma1Trace)
    .method("getPhiMixtureSigma2Trace", &Trace::getPhiMixtureSigma2Trace)
    .method("getPhiMixturePAcceptanceRateTrace",      &Trace::getPhiMixturePAcceptanceRateTrace)
    .method("getPhiMixtureMu1AcceptanceRateTrace",    &Trace::getPhiMixtureMu1AcceptanceRateTrace)
    .method("getPhiMixtureSigma1AcceptanceRateTrace", &Trace::getPhiMixtureSigma1AcceptanceRateTrace)
    .method("getPhiMixtureSigma2AcceptanceRateTrace", &Trace::getPhiMixtureSigma2AcceptanceRateTrace)

    // Phi-mixture-LN trace setters (task #12c.3): expose so
    // writeParameterObject / loadParameterObject can round-trip the
    // mixture-hyper traces through .Rdata serialisation.
    .method("setPhiMixturePTrace",      &Trace::setPhiMixturePTrace)
    .method("setPhiMixtureMu1Trace",    &Trace::setPhiMixtureMu1Trace)
    .method("setPhiMixtureSigma1Trace", &Trace::setPhiMixtureSigma1Trace)
    .method("setPhiMixtureSigma2Trace", &Trace::setPhiMixtureSigma2Trace)
    .method("setPhiMixturePAcceptanceRateTrace",      &Trace::setPhiMixturePAcceptanceRateTrace)
    .method("setPhiMixtureMu1AcceptanceRateTrace",    &Trace::setPhiMixtureMu1AcceptanceRateTrace)
    .method("setPhiMixtureSigma1AcceptanceRateTrace", &Trace::setPhiMixtureSigma1AcceptanceRateTrace)
    .method("setPhiMixtureSigma2AcceptanceRateTrace", &Trace::setPhiMixtureSigma2AcceptanceRateTrace)


    //Setter Functions:
    .method("setStdDevSynthesisRateTraces", &Trace::setStdDevSynthesisRateTraces)
    .method("setStdDevSynthesisRateAcceptanceRateTrace", &Trace::setStdDevSynthesisRateAcceptanceRateTrace)
    .method("setSynthesisRateTrace", &Trace::setSynthesisRateTrace)
    .method("setSynthesisRateAcceptanceRateTrace", &Trace::setSynthesisRateAcceptanceRateTrace)
    .method("setMixtureAssignmentTrace", &Trace::setMixtureAssignmentTrace)
    .method("setMixtureProbabilitiesTrace", &Trace::setMixtureProbabilitiesTrace)
    .method("setCodonSpecificAcceptanceRateTrace", &Trace::setCodonSpecificAcceptanceRateTrace)
    .method("setNseRateSpecificAcceptanceRateTrace", &Trace::setNseRateSpecificAcceptanceRateTrace)


    //ROC Specific:
    .method("getCodonSpecificParameterTraceByMixtureElementForCodon",
            &Trace::getCodonSpecificParameterTraceByMixtureElementForCodonR)
    .method("getSynthesisOffsetTrace", &Trace::getSynthesisOffsetTraceR)
    .method("getObservedSynthesisNoiseTrace", &Trace::getObservedSynthesisNoiseTraceR)
    .method("setSynthesisOffsetTrace", &Trace::setSynthesisOffsetTrace)
    .method("setSynthesisOffsetAcceptanceRateTrace", &Trace::setSynthesisOffsetAcceptanceRateTrace)
    .method("setObservedSynthesisNoiseTrace", &Trace::setObservedSynthesisNoiseTrace)
    .method("setCodonSpecificParameterTrace", &Trace::setCodonSpecificParameterTrace)

    //PANSE Specific
    .method("resizeNumberCodonSpecificParameterTrace", &Trace::resizeNumberCodonSpecificParameterTrace)
    .method("getPartitionFunctionTraces",&Trace::getPartitionFunctionTraces)
    .method("setPartitionFunctionTraces",&Trace::setPartitionFunctionTraces)


    //FONSE Specific
    .method("getInitiationCostTrace",&Trace::getInitiationCostTrace)
    .method("setInitiationCostTrace", &Trace::setInitiationCostTrace)
    ;
}
#endif

