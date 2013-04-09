import FWCore.ParameterSet.Config as cms

photonInvariantMassFilter = cms.EDFilter("PhotonInvariantMassFilter",
                                         photonCollection = cms.InputTag("goodPhotonsLowPtCut"),
                                         minMassCut = cms.double(45.),
                                         maxMassCut = cms.double(250.)
                                         )
