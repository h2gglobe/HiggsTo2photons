import FWCore.ParameterSet.Config as cms

invariantMassFilter = cms.EDFilter("InvariantMassFilter",
                                   electronCollection = cms.InputTag("goodElectronsOver5"),
                                   superClusterCollection = cms.InputTag("goodSCOver5"),
                                   minMassCut = cms.double(45),
                                   maxMassCut = cms.double(9999.)
                                   )
