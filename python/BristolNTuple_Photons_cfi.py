import FWCore.ParameterSet.Config as cms

rootTuplePhotons = cms.EDProducer("BristolNTuple_Photons",
#    InputTag = cms.InputTag('cleanPatPhotons'),
    InputTag = cms.InputTag('photons'),
    VertexTag = cms.InputTag('offlinePrimaryVertices'),
    ParticleFlowTag = cms.InputTag('particleFlow'),
    Prefix = cms.string('Photon'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(25),
)
