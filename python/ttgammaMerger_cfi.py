
import FWCore.ParameterSet.Config as cms

ttbarPhotonMerger2to5 = cms.EDFilter("TTGammaMerger",
    ptCut = cms.double(20.),
    drCut = cms.double(.1),
    filter = cms.bool(True),
)



