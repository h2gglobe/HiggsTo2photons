import FWCore.ParameterSet.Config as cms

process = cms.Process("H2gAnalyzer") 
process.load("Configuration.StandardSequences.Geometry_cff") 
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_35X_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# to run muon reco
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_3_5_4_RelValHGG_cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.p = cms.Path(process.h2ganalyzer)

process.h2ganalyzer.RootFileName = 'file:hgg.root'

process.h2ganalyzer.doReducedGen = False
process.h2ganalyzer.doGenParticles = False
process.h2ganalyzer.doSimHits = False
process.h2ganalyzer.doSimTracks = False
process.h2ganalyzer.doSimTrackPlusSimVertex = False
process.h2ganalyzer.doGenJet_mid = False
process.h2ganalyzer.doGenJet_it5 = False
process.h2ganalyzer.doGenJet_it7 = False

process.h2ganalyzer.doGenerator = False
process.h2ganalyzer.doReducedGen = False
process.h2ganalyzer.doGenParticles = False

process.h2ganalyzer.doGenJet_mid = False
process.h2ganalyzer.doGenJet_it5 = False
process.h2ganalyzer.doGenJet_it7 = False
    
process.h2ganalyzer.doSimHits = False
process.h2ganalyzer.doSimTracks = False
process.h2ganalyzer.doSimTrackPlusSimVertex = False
    
process.h2ganalyzer.doTrackingParticles = False

process.h2ganalyzer.doEcalRecHits = False

process.h2ganalyzer.doPreshowerHits = False
process.h2ganalyzer.doEcal = False
process.h2ganalyzer.doHcal = False
process.h2ganalyzer.doHFHcal = False
    
process.h2ganalyzer.doCaloTower = False
    
process.h2ganalyzer.doL1 = False
process.h2ganalyzer.doHLT = False
    
process.h2ganalyzer.doTracks = False
process.h2ganalyzer.doGsfTracks = False

process.h2ganalyzer.doVertices_std = False
process.h2ganalyzer.doVertices_pix = False
process.h2ganalyzer.doVtxCompat = False

process.h2ganalyzer.doElectron_std = False
process.h2ganalyzer.doMuon = False
process.h2ganalyzer.doPhoton = False
process.h2ganalyzer.doConvertedPhoton = False
process.h2ganalyzer.doLeptons = False

process.h2ganalyzer.doJet_mid = False
process.h2ganalyzer.doJet_it5 = False
process.h2ganalyzer.doJet_it7 = False
process.h2ganalyzer.doJet_it5pf = False 
process.h2ganalyzer.doJet_sis5pf = False
process.h2ganalyzer.doJet_kt4pf = False 

process.h2ganalyzer.doMet = False
process.h2ganalyzer.dotcMet = False
process.h2ganalyzer.doPFMet = False
process.h2ganalyzer.doHt = False

process.h2ganalyzer.doFastSim = False
process.h2ganalyzer.doAodSim  = False

#Minimal set

#process.h2ganalyzer.doL1 = True
#process.h2ganalyzer.doHLT = True
process.h2ganalyzer.doElectron_std = True
process.h2ganalyzer.doPhoton = True
process.h2ganalyzer.doConvertedPhoton = True
process.h2ganalyzer.doGsfTracks = True
process.h2ganalyzer.doPFMet = True
process.h2ganalyzer.doMet = True
process.h2ganalyzer.dotcMet = True

process.h2ganalyzer.Debug_Level = 0

process.GlobalTag.globaltag = 'MC_3XY_V25::All'

 
