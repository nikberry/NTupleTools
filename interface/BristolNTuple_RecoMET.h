#ifndef BristolNTupleRecoMETExtra
#define BristolNTupleRecoMETExtra

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class BristolNTuple_RecoMET : public edm::EDProducer {
 public:
  explicit BristolNTuple_RecoMET(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
};

#endif
