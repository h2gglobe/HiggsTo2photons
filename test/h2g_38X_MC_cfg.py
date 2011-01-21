# original version produced by job_maker on Mon Oct 18 04:40:41 2010
# command:
#   job_maker --data --skim2pho --name=Sept17ReReco --datasetpath=/EG/Run2010A-Sep17ReReco_v2/RECO --lumis=-1 --lumis_per_job=500 --outputdir=ssimon/data/2010/sept17rereco --json=~/json/goodrunlist_OCT15_json.txt


import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("H2gAnalyzer") 
process.load("Configuration.StandardSequences.GeometryDB_cff") 
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_38X_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process = cms.Process("H2gAnalyzer") 
#process.load("Configuration.StandardSequences.GeometryDB_cff") 
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_38X_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("PhysicsTools/PatAlgos/patSequences_cff")

#FLAGA (DATA TYPE)
flagData = 'OFF'
flagMC = 'ON'
flagFastSim = 'OFF'
flagNoSkim = 'OFF'

#check flags
flags = [flagData, flagMC, flagFastSim]
if (flags.count("ON") > 1):
    print "You are running too many things at the same time."
    sys.exit(-1)


#HLT label
hltLabel = "HLT"

#process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_3_8_3_RelValHGG_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200))
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
      '/store/mc/Fall10/GluGluToHToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/START38_V12-v2/0000/1C8E4A38-D8EB-DF11-924C-001F29079F98.root'
    )
)

#process.skipEvents = cms.untracked.PSet(input=cms.untracked.uint32(3500))
#skipEvents = cms.untracked.uint32(3500)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500



##################################### EVENT FILTER #######################################
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
process.L140or41 = process.hltLevel1GTSeed.clone()
process.L140or41.L1TechTriggerSeeding = cms.bool(True)
process.L140or41.L1SeedsLogicalExpression = cms.string('40 OR 41')

#process.hltPhysicsDeclared = cms.EDFilter("HLTHighLevel",
#                                          TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#                                          HLTPaths = cms.vstring("HLT_PhysicsDeclared"),
#                                          eventSetupPathsKey = cms.string(''),
#                                          andOr = cms.bool(True),
#                                          throw = cms.bool(True)
#                                          )

if (flagData == 'ON'):
    process.L10and34 = process.hltLevel1GTSeed.clone()
    process.L10and34.L1TechTriggerSeeding = cms.bool(True)
    process.L10and34.L1SeedsLogicalExpression = cms.string('0 AND 34 AND NOT (36,37,38,39)')
elif (flagMC == 'ON'):
    process.L10and34 = process.hltLevel1GTSeed.clone()
    process.L10and34.L1TechTriggerSeeding = cms.bool(True)
    process.L10and34.L1SeedsLogicalExpression = cms.string('34')


#process.hltJet15 = process.hltHighLevel.clone()
#process.hltJet15.TriggerResultsTag = cms.InputTag("TriggerResults","",hltLabel)
#process.hltJet15.HLTPaths = ("HLT_Jet15U",)
#process.hltJet15.andOr = True # True = OR, False = AND
#process.hltJet15.throw = False

#process.hltPhoton15 = process.hltHighLevel.clone()
#process.hltPhoton15.TriggerResultsTag = cms.InputTag("TriggerResults","",hltLabel)
#process.hltPhoton15.HLTPaths = ("HLT_Photon15_L1R", "HLT_Photon10_L1R", "HLT_Photon10_Cleaned_L1R", "HLT_Photon15_Cleaned_L1R", "HLT_Ele15_LW_L1R")
#process.hltPhoton15.andOr = True # True = OR, False = AND
#process.hltPhoton15.throw = False

#process.hltPhoton30 = process.hltHighLevel.clone()
#process.hltPhoton30.TriggerResultsTag = cms.InputTag("TriggerResults","",hltLabel)
#process.hltPhoton30.HLTPaths = ("HLT_Photon30_Cleaned_L1R", "HLT_Photon20_Cleaned_L1R")
#process.hltPhoton30.andOr = False # True = OR, False = AND
#process.hltPhoton30.throw = False

process.superClusterMerger =  cms.EDProducer("EgammaSuperClusterMerger",
                                             src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'),
                                                                 cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
                                             )

process.goodPhotonsLowPtCut = cms.EDFilter("PhotonSelector",
                                      src = cms.InputTag("photons"),
                                      cut = cms.string(
                                      "abs(superCluster.eta) < 2.5"
                                      " && superCluster.energy*sin(superCluster.position.theta) > 20."
                                      " && hadronicOverEm < 0.5 "
                                      " && trkSumPtHollowConeDR03 < 2.0*(3.5 + 0.001*superCluster.energy*sin(superCluster.position.theta))"
                                      " && ecalRecHitSumEtConeDR03 < 2.0*(4.2 + 0.006*superCluster.energy*sin(superCluster.position.theta))"
                                      " && hcalTowerSumEtConeDR03 < 2.0*(2.2 + 0.0025*superCluster.energy*sin(superCluster.position.theta))"
                                                      )
                                   )

process.TwoPhotonsLowPtCut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("goodPhotonsLowPtCut"),
                                       minNumber = cms.uint32(2)
                                     )

process.goodPhotonsHighPtCut = cms.EDFilter("PhotonSelector",
                                      src = cms.InputTag("photons"),
                                      cut = cms.string(
                                      "abs(superCluster.eta) < 2.5"
                                      " && superCluster.energy*sin(superCluster.position.theta) > 30."
                                      " && hadronicOverEm < 0.5 "
                                      " && trkSumPtHollowConeDR03 < 2.0*(3.5 + 0.001*superCluster.energy*sin(superCluster.position.theta))"
                                      " && ecalRecHitSumEtConeDR03 < 2.0*(4.2 + 0.006*superCluster.energy*sin(superCluster.position.theta))"
                                      " && hcalTowerSumEtConeDR03 < 2.0*(2.2 + 0.0025*superCluster.energy*sin(superCluster.position.theta))"
                                                      )
                                   )

process.OnePhotonsHighPtCut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("goodPhotonsHighPtCut"),
                                       minNumber = cms.uint32(1)
                                     )


process.noScraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15),
                                           maxd0 = cms.double(2)
                                           )

process.HFCoincidence =  cms.EDFilter("HFCoincidence",
                                          threshold = cms.untracked.double(2.0),
                                          HFRecHitCollection = cms.InputTag("hfreco"),
                                          #maskedChannels = cms.vint32( 8137, 8141, 8146, 8149, 8150, 8153 )
                                          )


if (flagData == 'ON'):
    process.goodEvents = cms.Sequence(process.noScraping * process.primaryVertexFilter)
    process.pathToCheck = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.HFCoincidence*process.L140or41)
    process.pathToCheck2 = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.L140or41*process.HFCoincidence)
elif (flagMC == 'ON'):
    process.goodEvents = cms.Sequence(process.noScraping * process.primaryVertexFilter)
    process.pathToCheck = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.HFCoincidence*process.L140or41)
    process.pathToCheck2 = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.L140or41*process.HFCoincidence)

if (flagNoSkim == 'ON'):
    process.eventFilter1 = cms.Sequence(process.dummySelector)
    process.eventFilter2 = cms.Sequence(process.dummySelector)
    process.eventFilter3 = cms.Sequence(process.dummySelector)
    process.eventFilter4 = cms.Sequence(process.dummySelector)
else:    
    process.eventFilter1 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut)
    process.eventFilter2 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut)
    process.eventFilter3 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut)
    process.eventFilter4 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut)
#process.eventFilter1 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut*process.goodPhotonsHighPtCut*process.OnePhotonsHighPtCut)
#process.eventFilter2 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut*process.goodPhotonsHighPtCut*process.OnePhotonsHighPtCut)
#process.eventFilter3 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut*process.goodPhotonsHighPtCut*process.OnePhotonsHighPtCut)
#process.eventFilter4 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut*process.goodPhotonsHighPtCut*process.OnePhotonsHighPtCut)
    

   
process.h2ganalyzer.RootFileName = 'file:hgg.root'

###################################DATA SEQUENCES#####################################
if(flagData == 'ON'):
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'])
   
    process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)

    process.p21 = cms.Path(process.eventFilter1*process.h2ganalyzerPath)
    process.p22 = cms.Path(process.eventFilter2*process.h2ganalyzerPath)
    process.p23 = cms.Path(process.eventFilter3*process.h2ganalyzerPath)
    process.p24 = cms.Path(process.eventFilter4*process.h2ganalyzerPath)


###################################MC SEQUENCES#####################################
if(flagMC == 'ON'):
    from PhysicsTools.PatAlgos.tools.coreTools import *
    
    process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)

    process.p21 = cms.Path(process.eventFilter1*process.h2ganalyzerPath)
    process.p22 = cms.Path(process.eventFilter2*process.h2ganalyzerPath)
    process.p23 = cms.Path(process.eventFilter3*process.h2ganalyzerPath)
    process.p24 = cms.Path(process.eventFilter4*process.h2ganalyzerPath)
    
        
###################################FASTSIM SEQUENCES#####################################
if(flagFastSim == 'ON'):
    process.simulation = cms.Sequence(process.ProductionFilterSequence*process.simulationWithFamos)
    process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)
    process.p2 = cms.Path(process.simulation+process.reconstructionWithFamos+process.eventFilter1*process.h2ganalyzerPath)


###################################H2GANALYZER PARAMETERS#####################################

# RUNNING MODULES
process.h2ganalyzer.doGenerator                 = False
process.h2ganalyzer.doGenParticles              = True
process.h2ganalyzer.doGenJet_mid                = True
process.h2ganalyzer.doGenJet_it5                = True
process.h2ganalyzer.doGenJet_it7                = True
process.h2ganalyzer.doReducedGen                = False
process.h2ganalyzer.doSimHits                   = False
process.h2ganalyzer.doSimTracks                 = True
process.h2ganalyzer.doSimTrackPlusSimVertex     = True
process.h2ganalyzer.doTrackingParticles         = False
process.h2ganalyzer.doEcalRecHits               = True
process.h2ganalyzer.doPreshowerHits             = True
process.h2ganalyzer.doEcal                      = True
process.h2ganalyzer.doHcal                      = True
process.h2ganalyzer.doHFHcal                    = True
process.h2ganalyzer.doCaloTower                 = True
process.h2ganalyzer.doL1                        = False
process.h2ganalyzer.doHLT                       = True
process.h2ganalyzer.doTracks                    = True
process.h2ganalyzer.doGsfTracks                 = True
process.h2ganalyzer.doVertices_std              = True
process.h2ganalyzer.doVertices_pix              = True
process.h2ganalyzer.doVtxCompat                 = False
process.h2ganalyzer.doElectron_std              = True
process.h2ganalyzer.doMuon                      = True
process.h2ganalyzer.doPhoton                    = True
process.h2ganalyzer.doConvertedPhoton           = False
process.h2ganalyzer.doLeptons                   = True
process.h2ganalyzer.doJet_mid                   = True
process.h2ganalyzer.doJet_it5                   = True
process.h2ganalyzer.doJet_it7                   = True
process.h2ganalyzer.doJet_it5pf                 = True
process.h2ganalyzer.doJet_sis5pf                = True
process.h2ganalyzer.doJet_kt4pf                 = True
process.h2ganalyzer.doMet                       = True
process.h2ganalyzer.dotcMet                     = True
process.h2ganalyzer.doPFMet                     = True
process.h2ganalyzer.doHt                        = True



if(flagData == 'ON'):
    process.h2ganalyzer.doGenerator                 = False
    process.h2ganalyzer.doReducedGen                = False
    process.h2ganalyzer.doGenParticles              = False
    process.h2ganalyzer.doGenJet_mid                = False
    process.h2ganalyzer.doGenJet_it5                = False
    process.h2ganalyzer.doGenJet_it7                = False
    process.h2ganalyzer.doSimHits                   = False
    process.h2ganalyzer.doSimTracks                 = False
    process.h2ganalyzer.doSimTrackPlusSimVertex     = False
    process.h2ganalyzer.doTrackingParticles         = False
    process.h2ganalyzer.doL1 = False
    process.h2ganalyzer.doHLT = True
    process.h2ganalyzer.Debug_Level = 0
    process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag("TriggerResults","",hltLabel)
    process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
    process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","",hltLabel)
    process.GlobalTag.globaltag = 'GR_R_38X_V12::All' #CHECK ALWAYS
    
elif(flagMC == 'ON'):
    process.h2ganalyzer.doSimHits = False
    process.h2ganalyzer.doTrackingParticles = False
    process.h2ganalyzer.doL1 = False
    process.h2ganalyzer.doHLT = True
    process.GlobalTag.globaltag = 'START38_V13::All' #CHECK ALWAYS
    process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag("TriggerResults","",hltLabel)
    process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
    process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","",hltLabel)

elif(flagFastSim == 'ON'):
    process.h2ganalyzer.doSimHits = False
    process.h2ganalyzer.doSimTracks = False
    process.h2ganalyzer.doSimTrackPlusSimVertex = False
    process.h2ganalyzer.doTrackingParticles = False
    process.h2ganalyzer.doL1 = False
    process.h2ganalyzer.doHLT = False
    process.h2ganalyzer.doJet_mid = False
    process.h2ganalyzer.doJet_it5 = False
    process.h2ganalyzer.doJet_it7 = False
    process.h2ganalyzer.doFastSim = True
    process.GlobalTag.globaltag = "STARTUP38X_V13::All"  #CHECK ALWAYS



