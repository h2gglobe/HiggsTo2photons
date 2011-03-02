import FWCore.ParameterSet.Config as cms

process = cms.Process("Globe") 
process.load("Configuration.StandardSequences.GeometryDB_cff") 
process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_39X_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

#from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *

process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_RelValDUMMY_cfi")
#process.skipEvents = cms.untracked.PSet(input=cms.untracked.uint32(3500))
#skipEvents = cms.untracked.uint32(3500)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

# if wSelect is loaded after wFilter there is an error.
#FIXME ADD PHOTON SKIM
#process.goodEvents = cms.Sequence(process.noScraping * process.primaryVertexFilter)
#process.pathToCheck = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.HFCoincidence*process.L140or41)
#process.pathToCheck2 = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.L140or41*process.HFCoincidence)

process.h2ganalyzer.RootFileName = 'reco_data_test.root'
process.h2ganalyzer.Debug_Level = 0

process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)

process.p11 = cms.Path(process.h2ganalyzerPath)

process.h2ganalyzer.doGenJet_algo1 = False
process.h2ganalyzer.doGenJet_algo2 = False
process.h2ganalyzer.doGenJet_algo3 = False
process.h2ganalyzer.doGenParticles = False
process.h2ganalyzer.doReducedGen = False

process.h2ganalyzer.doSimTrackPlusSimVertex = False

process.h2ganalyzer.doL1 = True
process.h2ganalyzer.doHLT = True
process.GlobalTag.globaltag = "START39_v8::All"
process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag("TriggerResults","", "HLT")
process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","", "HLT")
