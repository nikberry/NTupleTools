#ifndef BristolNTuplePhotons
#define BristolNTuplePhotons

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"

class BristolNTuple_Photons: public edm::EDProducer {
public:
    explicit BristolNTuple_Photons(const edm::ParameterSet&);

private:
    void produce(edm::Event &, const edm::EventSetup &);
    const edm::InputTag inputTag;
    const edm::InputTag vertexTag;
    const edm::InputTag particleFlowTag;
    const std::string prefix, suffix;
    const unsigned int maxSize;
    PFIsolationEstimator isolator;
};

#endif
