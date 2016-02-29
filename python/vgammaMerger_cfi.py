
import FWCore.ParameterSet.Config as cms

vplusphotonMerger2to7 = cms.EDFilter("VGammaMerger",                                                                                                                                                               
    ptCut = cms.double(13.),
    drCut = cms.double(.3),
    etaCut = cms.double(3),
    filter = cms.bool(True),
)

