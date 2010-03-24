import FWCore.ParameterSet.Config as cms

process = cms.Process("Globe") 
process.load("Configuration.StandardSequences.Geometry_cff") 
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.load("HiggsAnalysis.HiggsTo2photons.globeanalyzer_33X_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# to run muon reco
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_3_1_2_RelValHGG_cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.p = cms.Path(process.globe)

process.globe.RootFileName = 'file:hgg.root'
process.globe.doTrackingParticles = False
process.globe.doSimHits = False
process.globe.doSimTracks = True
process.globe.doSimTrackPlusSimVertex = True
process.globe.doJet_it5pf = True
process.globe.doJet_sis5pf = True
process.globe.doJet_kt4pf = True
process.globe.dotcMet = True
process.globe.doPFMet = True
process.globe.doL1 = False
process.globe.doHLT = True
process.globe.Debug_Level = 0

process.GlobalTag.globaltag = 'STARTUP31X_V2::All'

 
