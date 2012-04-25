import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *

hggPhotonIDConfiguration = cms.PSet(hggPhotonIDCuts)

cicFilter = cms.EDFilter("cicFilter",
                         nPhotons = cms.int32(1),
                         PhotonCollection = cms.InputTag("photons"),
                         PFCollection = cms.InputTag("particleFlow"),
                         ElectronCollection = cms.InputTag("gsfElectrons"),
                         VertexCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
                         TrackCollection = cms.InputTag("generalTracks"),
                         RhoCollection = cms.InputTag("kt6PFJets","rho"),
                         CutLevel = cms.int32(11), # 11 means preselection
                         useOR = cms.bool(True),
                         hggPhotonIDConfiguration = cms.PSet(hggPhotonIDCuts)
                         )


cicFilterSequence = cms.Sequence(cicFilter)
