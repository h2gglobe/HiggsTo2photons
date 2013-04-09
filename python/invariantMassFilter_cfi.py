import FWCore.ParameterSet.Config as cms

invariantMassFilter = cms.EDFilter("InvariantMassFilter",
                                   electronCollection = cms.InputTag("goodElectronsOver5"),
                                   superClusterCollection = cms.InputTag("goodSCOver5"),
                                   ElectronElectronInvMass = cms.bool(False),
                                   minMassCut = cms.double(50),
                                   maxMassCut = cms.double(9999.)
                                   )

