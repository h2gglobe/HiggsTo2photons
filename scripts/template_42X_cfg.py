import FWCore.ParameterSet.Config as cms
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *

#DATA TYPE
flagData = 'OFF'
flagMC = 'OFF'
flagFastSim = 'ON'

#SKIM TYPE
flagSkimDiphoton = 'OFF'
flagVLPreSelection = 'OFF'
flagNoSkim = 'OFF'
flagMMgSkim = 'OFF'
flagSkimworz = 'OFF'
flagSkim1El = 'OFF'
flagAddPdfWeight = 'OFF'

#ADDITIONAL OPTIONS
flagAOD = 'ON'
jobMaker = 'jobmaker unknown'

if (not((flagNoSkim is 'ON') ^ (flagSkimDiphoton is 'ON') ^ (flagMMgSkim is 'ON') ^ (flagVLPreSelection is 'ON') ^ (flagSkim1El is 'ON') ^ (flagSkimworz is 'ON'))):
  print "You must skim or not skim... these are your options"
  exit(-1)

process = cms.Process("Globe") 
process.load("Configuration.StandardSequences.GeometryDB_cff") 
process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_42X_cfi")
#pi0 disc
process.load("RecoEcal.EgammaClusterProducers.preshowerClusterShape_cfi")
process.load("EgammaAnalysis.PhotonIDProducers.piZeroDiscriminators_cfi")
  
if flagAOD is 'OFF':
  #rerun ConvId
  process.load("RecoEgamma.EgammaPhotonProducers.conversionTrackSequence_cff")
  # NOTICE: You need the following two python files to rerun the conversion tracking with ECAL association
  process.load("RecoEgamma.EgammaPhotonProducers.conversionTrackCandidates_cfi")
  process.load("RecoEgamma.EgammaPhotonProducers.ckfOutInTracksFromConversions_cfi")
  ###############

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Mixing')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('HiggsAnalysis.HiggsTo2photons.ZMuSkim_cff')
process.load('HiggsAnalysis.HiggsTo2photons.photonReRecoForMMG_cfi')

if flagSkimDiphoton == 'ON':
#  process.load('Configuration.Skimming.PDWG_DiPhoton_SD_cff')
     process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
     process.DiPhotonHltFilter = copy.deepcopy(hltHighLevel)
     process.DiPhotonHltFilter.throw = cms.bool(False)
     process.DiPhotonHltFilter.HLTPaths = ["HLT_Photon*_CaloId*_Iso*_Photon*_CaloId*_Iso*_*","HLT_Photon*_R9Id*_Photon*_R9Id*_*","HLT_Photon*_R9Id*_Photon*_CaloId*_Iso*_*","HLT_Photon*_CaloId*_Iso*_Photon*_R9Id*_*"]



process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_RelValDUMMY_cfi")
#process.skipEvents = cms.untracked.PSet(input=cms.untracked.uint32(3500))
#skipEvents = cms.untracked.uint32(3500)

hltLabel = "HLT"

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500


process.superClusterMerger =  cms.EDProducer("EgammaSuperClusterMerger",
                                             src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'),
                                                                 cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
                                             )


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15),
                                           maxd0 = cms.double(2)
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


process.dummySelector = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("gsfElectrons"),
                                     minNumber = cms.uint32(0)
                                     )

process.goodElectronsOver5 = cms.EDFilter("GsfElectronSelector",
                                          filter = cms.bool(True),
                                          src = cms.InputTag("gsfElectrons"),
                                          cut = cms.string('superCluster().rawEnergy()*sin((2*atan(exp(superCluster().eta())))) > 5.')
                                          )

process.superClusterMerger =  cms.EDProducer("EgammaSuperClusterMerger",
                                             src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'),
                                                                 cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
                                             )

process.goodSCOver5 = cms.EDFilter("SuperClusterSelector",
                                   filter = cms.bool(True),
                                   src = cms.InputTag("superClusterMerger"),
                                   cut = cms.string('rawEnergy()*sin((2*atan(exp(eta())))) > 5.')
                                   )

process.pdfWeights = cms.EDProducer("PdfWeightProducer",
                                    #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
                                    #GenTag = cms.untracked.InputTag("genParticles"),
                                    PdfInfoTag = cms.untracked.InputTag("generator"),
                                    PdfSetNames = cms.untracked.vstring("cteq66.LHgrid")
                                    )

#process.goodEvents = cms.Sequence(process.noScraping * process.primaryVertexFilter)
#process.pathToCheck = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.HFCoincidence*process.L140or41)
#process.pathToCheck2 = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.L140or41*process.HFCoincidence)

process.diMuonSelSeq.remove(process.ZMuHLTFilter)

if flagVLPreSelection == 'ON':
  process.eventFilter1 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut) # for bkg
  process.eventFilter2 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut) # for bkg
elif flagSkimDiphoton == 'ON':
  process.eventFilter1 = cms.Sequence(process.DiPhotonHltFilter) # for some data
  process.eventFilter2 = cms.Sequence(process.DiPhotonHltFilter)      # for some data
elif flagNoSkim == 'ON':    
  process.eventFilter1 = cms.Sequence(process.dummySelector)   #for signal MC
  process.eventFilter2 = cms.Sequence(process.dummySelector)   #for signal MC
elif flagMMgSkim == 'ON':
  process.eventFilter1 = cms.Sequence(process.diMuonSelSeq*process.photonReReco)
  process.eventFilter2 = cms.Sequence(process.diMuonSelSeq*process.photonReReco)
elif flagSkimworz == 'ON':
  process.load('HiggsAnalysis.HiggsTo2photons.eidFilter_cfi')
  process.load('HiggsAnalysis.HiggsTo2photons.invariantMassFilter_cfi')
  process.eventFilter1 = cms.Sequence(process.goodElectronsOver5*process.superClusterMerger*process.goodSCOver5*process.invariantMassFilter*process.electronIdentificationFilter)
  process.eventFilter2= cms.Sequence(process.goodElectronsOver5*process.superClusterMerger*process.goodSCOver5*process.invariantMassFilter*process.electronIdentificationFilter)
elif flagSkim1El == 'ON':
  process.eventFilter1 = cms.Sequence(process.goodElectronsOver5)
  process.eventFilter2 = cms.Sequence(process.goodElectronsOver5)


process.h2ganalyzer.RootFileName = 'aod_mc_test.root'
process.h2ganalyzer.Debug_Level = 0

##-------------------- PFNoPU for PF Isolation Electrons -------------
process.load("CommonTools.ParticleFlow.pfPileUp_cfi")
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets = process.kt6PFJets.clone(rParam = 0.6, doRhoFastjet = True)
process.ak5PFJets.doAreaFastjet = True
process.kt6PFJetsForRhoCorrection = process.kt6PFJets.clone(rParam = 0.6, doRhoFastjet = True)
process.kt6PFJetsForRhoCorrection.Rho_EtaMax = cms.double(2.5)


##-------------------- Filter to skip bugged events with non conserved energy -------
process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")

# event counters
process.processedEvents = cms.EDProducer("EventCountProducer")
process.eventCounters = cms.Sequence(process.processedEvents)

if (flagFastSim == 'OFF'):
  process.eventCounters = cms.Sequence(process.totalKinematicsFilter * process.processedEvents)
  
if (flagAddPdfWeight == 'ON'):
  process.eventCounters *= cms.Sequence(process.pdfWeights)
  
process.h2ganalyzer.globalCounters.extend(['processedEvents']) 

# PFIsolation photons
process.load("HiggsAnalysis.HiggsTo2photons.pfIsolation_cff")

process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)

#################################################
# Define path, first for AOD case then for RECO #
#################################################
#process.pfBasedPhotonIsoSequence
#process.pfSelectedPhotons

process.p11 = cms.Path(process.eventCounters*process.eventFilter1*process.pfPileUp)
  
if (flagFastSim == 'OFF'):
  process.p11 *= process.piZeroDiscriminators
    
process.p11 *= (process.kt6PFJets* process.ak5PFJets* process.kt6PFJetsForRhoCorrection* process.h2ganalyzerPath)
  
process.p12 = copy.deepcopy(process.p11)
process.p12.replace(process.eventFilter1, process.eventFilter2)

if (flagAOD is 'OFF'):
  process.p11.insert(-1, (process.conversionTrackCandidates*process.ckfOutInTracksFromConversions*process.preshowerClusterShape*process.piZeroDiscriminators))

  process.p12.insert(-1, (process.conversionTrackCandidates*process.ckfOutInTracksFromConversions*process.preshowerClusterShape*process.piZeroDiscriminators))

#################################################
# End of Path definition                        #
#################################################

process.h2ganalyzer.JobMaker = jobMaker

if (flagAddPdfWeight == 'ON'):
  process.h2ganalyzer.doPdfWeight = True 

if (flagFastSim is 'ON'):
  process.h2ganalyzer.doFastSim = True

if (flagMC is 'ON' and flagFastSim is 'ON'):
  process.h2ganalyzer.doGenJet_algo1 = False
  process.h2ganalyzer.doGenJet_algo2 = False
  process.h2ganalyzer.doGenJet_algo3 = False
  process.h2ganalyzer.doGenParticles = False
  process.h2ganalyzer.doGenMet = False
  process.h2ganalyzer.doReducedGen = False
  process.h2ganalyzer.doGenVertices = False
elif (flagMC is 'ON' and flagFastSim is 'OFF'):
  process.h2ganalyzer.doGenJet_algo1 = True
  process.h2ganalyzer.doGenJet_algo2 = True
  process.h2ganalyzer.doGenJet_algo3 = True
  process.h2ganalyzer.doGenParticles = True
  process.h2ganalyzer.doGenMet = True
  process.h2ganalyzer.doReducedGen = True
  process.h2ganalyzer.doGenVertices = True
elif flagData is 'ON':
  process.h2ganalyzer.doPileup = False
  process.h2ganalyzer.doGenJet_algo1 = False
  process.h2ganalyzer.doGenJet_algo2 = False
  process.h2ganalyzer.doGenJet_algo3 = False
  process.h2ganalyzer.doGenParticles = False
  process.h2ganalyzer.doGenVertices = False
  process.h2ganalyzer.doGenMet = False
  process.h2ganalyzer.doReducedGen = False

if flagMC is 'ON' and flagAOD is 'OFF':
  process.h2ganalyzer.doSimTracks = True
  process.h2ganalyzer.doSimTrackPlusSimVertex = False

if flagAOD is 'ON':
  process.h2ganalyzer.doAodSim = True
  process.h2ganalyzer.doHcal = False
  process.h2ganalyzer.doHFHcal = False
  process.h2ganalyzer.doPreshowerHits = False
else:
  process.h2ganalyzer.doAodSim = False
  process.h2ganalyzer.doHcal = True
  process.h2ganalyzer.doHFHcal = True
  process.h2ganalyzer.doPreshowerHits = True
  process.h2ganalyzer.EcalHitEBColl = cms.InputTag("ecalRecHit","EcalRecHitsEB")
  process.h2ganalyzer.EcalHitEEColl = cms.InputTag("ecalRecHit","EcalRecHitsEE")
  process.h2ganalyzer.HcalHitsBEColl = cms.InputTag("hbhereco")
  process.h2ganalyzer.HcalHitsFColl = cms.InputTag("hfreco")
  process.h2ganalyzer.HcalHitsHoColl = cms.InputTag("horeco")
  process.h2ganalyzer.BarrelBasicClusterColl = cms.InputTag("")
  process.h2ganalyzer.BarrelBasicClusterShapeColl = cms.InputTag("multi5x5BasicClusters","multi5x5BarrelShapeAssoc")
  process.h2ganalyzer.JetTrackAssociationColl_algo3 = cms.InputTag("kt4JetTracksAssociatorAtVertex")


process.h2ganalyzer.doL1 = True
process.h2ganalyzer.doHLT = True

process.GlobalTag.globaltag = "START39_v8::All"
process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag("TriggerResults","", hltLabel)
process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","", hltLabel)
