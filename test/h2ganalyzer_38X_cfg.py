import FWCore.ParameterSet.Config as cms

process = cms.Process("H2gAnalyzer") 
process.load("Configuration.StandardSequences.GeometryDB_cff") 
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_38X_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

#process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_3_8_3_RelValHGG_cfi")
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
      '/store/mc/Fall10/GluGluToHToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/START38_V12-v2/0000/1C8E4A38-D8EB-DF11-924C-001F29079F98.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.p = cms.Path(process.h2ganalyzer)

process.h2ganalyzer.RootFileName = 'file:hgg.root'
process.h2ganalyzer.doTrackingParticles = False
process.h2ganalyzer.doSimHits = False
process.h2ganalyzer.doSimTracks = True
process.h2ganalyzer.doSimTrackPlusSimVertex = True
process.h2ganalyzer.doJet_it5pf = True
process.h2ganalyzer.doJet_sis5pf = True
process.h2ganalyzer.doJet_kt4pf = True
process.h2ganalyzer.dotcMet = True
process.h2ganalyzer.doPFMet = True
process.h2ganalyzer.doL1 = False
process.h2ganalyzer.doHLT = True
process.h2ganalyzer.Debug_Level = 0

process.GlobalTag.globaltag = 'START38_V9::All'

 
