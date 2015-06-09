
import FWCore.ParameterSet.Config as cms

ttbarPhotonMerger2to5 = cms.EDFilter("TTGammaMerger",
    ptCut = cms.double(13.),
    drCut = cms.double(.3),
    etaCut = cms.double(3),
    filter = cms.bool(True),
)



