import FWCore.ParameterSet.Config as cms
import copy

#DATA TYPE
flagData = 'OFF'
flagMC = 'OFF'
flagFastSim = 'OFF'
flagPG = 'OFF'

#SKIM TYPE
flagSkimDiphoton = 'OFF'
flagSkimPJet = 'OFF'
flagVLPreSelection = 'OFF'
flagMyPreSelection = 'OFF'
flagNoSkim = 'OFF'
flagMuMuSkim = 'OFF'
flagMMgSkim = 'OFF'
flagSkimworz = 'OFF'
flagSkim1El = 'OFF'
flagAddPdfWeight = 'OFF'
flagSkimHmm = 'OFF'
flagSkimHee = 'OFF'
flagSkimMu = 'OFF'

#ADDITIONAL OPTIONS
flagAOD = 'ON'
jobMaker = 'jobmaker unknown'

if (not((flagNoSkim is 'ON') ^ (flagSkimDiphoton is 'ON') ^ (flagMMgSkim is 'ON') ^ (flagVLPreSelection is 'ON') ^ (flagSkim1El is 'ON') ^ (flagSkimworz is 'ON') ^ (flagMyPreSelection is 'ON') ^ (flagSkimPJet is 'ON')^ (flagSkimHmm is 'ON')^ (flagSkimHee is 'ON') ^ (flagSkimMu is 'ON'))):
  print "You must skim or not skim... these are your options"
  exit(-1)

process = cms.Process("Globe")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_53X_cfi")
#pi0 disc
process.load("RecoEcal.EgammaClusterProducers.preshowerClusterShape_cfi")
process.load("EgammaAnalysis.PhotonIDProducers.piZeroDiscriminators_cfi")
#FIXME
process.piZeroDiscriminators.preshClusterShapeProducer = cms.string('multi5x5PreshowerClusterShape')
process.piZeroDiscriminators.preshClusterShapeCollectionX = cms.string('multi5x5PreshowerXClustersShape')
process.piZeroDiscriminators.preshClusterShapeCollectionY = cms.string('multi5x5PreshowerYClustersShape')

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
process.load('HiggsAnalysis.HiggsTo2photons.cicFilter_cfi')

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedElectrons = cms.PSet(
  initialSeed = cms.untracked.uint32(1),
  engineName = cms.untracked.string('TRandom3')
  ),
                                                   )

if flagSkimMu == 'ON':
  process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
  process.HLTSingleMu = copy.deepcopy(process.hltHighLevel)
  process.HLTSingleMu.throw = cms.bool(False)
  process.HLTSingleMu.HLTPaths = ["HLT_IsoMu24_eta2p1_*",]

if flagSkimworz == 'ON':
  process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
  process.TPHltFilter = copy.deepcopy(process.hltHighLevel)
  process.TPHltFilter.throw = cms.bool(False)
  process.TPHltFilter.HLTPaths = ["HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_*",
                                  "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_*",
                                  "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_*"]

if flagSkimDiphoton == 'ON':
  process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
  process.DiPhotonHltFilter = copy.deepcopy(process.hltHighLevel)
  process.DiPhotonHltFilter.throw = cms.bool(False)
  process.DiPhotonHltFilter.HLTPaths = ["HLT_Photon*_CaloId*_Iso*_Photon*_CaloId*_Iso*_*","HLT_Photon*_CaloId*_Iso*_Photon*_R9Id*_*","HLT_Photon*_R9Id*_Photon*_CaloId*_Iso*_*","HLT_Photon*_R9Id*_Photon*_R9Id*_*","HLT_Photon*_R9Id*_OR_CaloId*_Iso*_Photon*_R9Id*_OR_CaloId*_Iso*_*","HLT_Photon*_R9Id*_OR_CaloId*_Iso*_Photon*_*"]
  #process.load('Configuration.Skimming.PDWG_DiPhoton_SD_cff')

if flagSkimPJet == 'ON':
  process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
  process.PhotonHltFilter = copy.deepcopy(process.hltHighLevel)
  process.PhotonHltFilter.throw = cms.bool(False)
  process.PhotonHltFilter.HLTPaths = ["HLT_Photon*_CaloIdVL_IsoL_v*"]

process.load("HiggsAnalysis.HiggsTo2photons.photonInvariantMassFilter_cfi")
process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_RelValDUMMY_cfi")
#process.source.skipEvents = cms.untracked.uint32(3500)

hltLabel = "HLT"

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  SkipEvent = cms.untracked.vstring('FatalRootError','InvalidReference')
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1001))

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

process.goodSuperclusters = cms.EDFilter("PhotonSelector",
                                         src = cms.InputTag("photons"),
                                         cut = cms.string(
                                         "abs(superCluster.eta) < 2.5"
                                         " && superCluster.energy*sin(superCluster.position.theta) > 20."
                                         )
                                         )

process.OnePhotonsHighPtCut = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("goodPhotonsHighPtCut"),
                                           minNumber = cms.uint32(1)
                                           )

process.OneSupercluster = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("goodSuperclusters"),
                                       minNumber = cms.uint32(1)
                                       )


process.dummySelector = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("gsfElectrons"),
                                     minNumber = cms.uint32(0)
                                     )

process.muonsOver20 = cms.EDFilter("MuonSelector",
                                   filter = cms.bool(True),
                                   src = cms.InputTag("muons"),
                                   cut = cms.string('pt() > 20.')
                                   )

process.muonCounter = cms.EDFilter("CandViewCountFilter",
                                   src = cms.InputTag("muonsOver20"),
                                   minNumber = cms.uint32(2)
                                   )

process.electronsOver20 = cms.EDFilter("GsfElectronSelector",
                                       filter = cms.bool(True),
                                       src = cms.InputTag("gsfElectrons"),
                                       cut = cms.string('superCluster().rawEnergy()*sin((2*atan(exp(superCluster().eta())))) > 20.')
                                       )

process.electronCounter = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("electronsOver20"),
                                       minNumber = cms.uint32(2)
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

process.eventFilter1 = cms.Sequence()
process.eventFilter2 = cms.Sequence()

if flagVLPreSelection == 'ON':
  process.eventFilter1 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut+process.photonInvariantMassFilter) # for bkg
  process.eventFilter2 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut+process.photonInvariantMassFilter) # for bkg
elif flagMyPreSelection == 'ON':
  process.eventFilter1 = cms.Sequence(process.cicFilterSequence)
  process.eventFilter2 = cms.Sequence(process.cicFilterSequence)
elif flagSkimDiphoton == 'ON':
  process.eventFilter1 = cms.Sequence(process.DiPhotonHltFilter) # for some data
  process.eventFilter2 = cms.Sequence(process.DiPhotonHltFilter)      # for some data
  #process.eventFilter1 = cms.Sequence(process.CaloIdIsoPhotonPairsFilter) # for some data
  #process.eventFilter2 = cms.Sequence(process.R9IdPhotonPairsFilter)      # for some data
elif flagData == 'ON' and flagSkimPJet == 'ON':
  process.eventFilter1 = cms.Sequence(process.PhotonHltFilter)
  process.eventFilter2 = cms.Sequence(process.PhotonHltFilter)
elif flagMC == 'ON' and flagSkimPJet == 'ON':
  process.eventFilter1 = cms.Sequence(process.PhotonHltFilter*process.superClusterMerger*process.goodPhotonsHighPtCut*process.OnePhotonsHighPtCut) # for bkg
  process.eventFilter2 = cms.Sequence(process.PhotonHltFilter*process.superClusterMerger*process.goodPhotonsHighPtCut+process.OnePhotonsHighPtCut) # for bkg
elif flagNoSkim == 'ON':
  process.eventFilter1 = cms.Sequence(process.dummySelector)   #for signal MC
  process.eventFilter2 = cms.Sequence(process.dummySelector)   #for signal MC
elif flagMuMuSkim == 'ON':
  process.dimuons.cut = 'mass > 60'
  process.eventFilter1 = cms.Sequence(process.diMuonSelSeq)
  process.eventFilter2 = cms.Sequence(process.diMuonSelSeq)
elif flagMMgSkim == 'ON':
  process.eventFilter1 = cms.Sequence(process.diMuonSelSeq*process.photonReReco)
  process.eventFilter2 = cms.Sequence(process.diMuonSelSeq*process.photonReReco)
elif flagSkimHmm == 'ON':
  process.eventFilter1 = cms.Sequence(process.muonsOver20*process.muonCounter)
  process.eventFilter2 = cms.Sequence(process.muonsOver20*process.muonCounter)
elif flagSkimHee == 'ON':
  process.eventFilter1 = cms.Sequence(process.electronsOver20*process.electronCounter)
  process.eventFilter2 = cms.Sequence(process.electronsOver20*process.electronCounter)
elif flagSkimworz == 'ON':
  if (flagData == 'ON'):
    process.eventFilter1 = cms.Sequence(process.TPHltFilter)
    process.eventFilter2 = cms.Sequence(process.TPHltFilter)
  else:
    process.eventFilter1 = cms.Sequence()
    process.eventFilter2 = cms.Sequence()
  #process.load('HiggsAnalysis.HiggsTo2photons.eidFilter_cfi')
  process.load('HiggsAnalysis.HiggsTo2photons.invariantMassFilter_cfi')
  process.eventFilter1 *= cms.Sequence(process.goodElectronsOver5*process.superClusterMerger*process.goodSCOver5*process.invariantMassFilter)
  process.eventFilter2 *= cms.Sequence(process.goodElectronsOver5*process.superClusterMerger*process.goodSCOver5*process.invariantMassFilter)
elif flagSkim1El == 'ON':
  process.eventFilter1 = cms.Sequence(process.goodElectronsOver5)
  process.eventFilter2 = cms.Sequence(process.goodElectronsOver5)
elif flagSkimMu == 'ON':
  process.eventFilter1 = cms.Sequence(process.HLTSingleMu*process.superClusterMerger*process.goodSuperclusters*process.OneSupercluster)
  process.eventFilter2 = cms.Sequence(process.HLTSingleMu*process.superClusterMerger*process.goodSuperclusters*process.OneSupercluster)


process.h2ganalyzer.RootFileName = 'aod_mc_test.root'
process.h2ganalyzer.Debug_Level = 0

##---------------------ELECTRON REGRESSION AND SMEARING ------------------------------
process.load("EgammaAnalysis.ElectronTools.calibratedElectrons_cfi")

# dataset to correct
if (flagMC == 'ON'):
  process.calibratedElectrons.isMC = cms.bool(True)
  process.calibratedElectrons.inputDataset = cms.string("Summer12_DR53X_HCP2012")
else:
  process.calibratedElectrons.isMC = cms.bool(False)
  process.calibratedElectrons.inputDataset = cms.string("Moriond2013")

process.calibratedElectrons.updateEnergyError = cms.bool(True)
process.calibratedElectrons.applyCorrections = cms.int32(1)
process.calibratedElectrons.smearingRatio = cms.double(0.607)
process.calibratedElectrons.verbose = cms.bool(False)
#process.calibratedElectrons.synchronization = cms.bool(True)

process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('gsfElectrons')
process.eleRegressionEnergy.inputCollectionType = cms.uint32(0)
process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
process.eleRegressionEnergy.produceValueMaps = cms.bool(True)

##-------------------- ANOMALOUS HCAL LASER CORRECTION FILTER ------------------------
#process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cff")
##-------------------- ANOMALOUS ECAL LASER CORRECTION FILTER ------------------------
#process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
##-------------------- PFIsolation for Electrons -------------------------------------
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
#process.phoIsoSequence = setupPFPhotonIso(process, 'photons')
##-------------------- PFNoPU for PF Isolation Electrons -----------------------------
process.load("CommonTools.ParticleFlow.pfPileUp_cfi")
process.pfPileUp.PFCandidates = cms.InputTag("particleFlow")
##-------------------- Import the JEC services ---------------------------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------------------
#process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets = process.kt6PFJets.clone(rParam = 0.6, doRhoFastjet = True)
#process.ak5PFJets.doAreaFastjet = True
#process.kt6PFJetsForRhoCorrection = process.kt6PFJets.clone(rParam = 0.6, doRhoFastjet = True)
#process.kt6PFJetsForRhoCorrection.Rho_EtaMax = cms.double(2.5)

##-------------------- Filter to skip bugged events with non conserved energy --------
process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")

# event counters
process.processedEvents = cms.EDProducer("EventCountProducer")

process.eventCounters = cms.Sequence(process.processedEvents)

if (flagFastSim == 'OFF' and flagMC == 'ON' and flagPG == 'OFF'):
  process.eventCounters = cms.Sequence(process.totalKinematicsFilter * process.processedEvents)

if (flagAddPdfWeight == 'ON' and flagMC == 'ON'):
  process.eventCounters *= cms.Sequence(process.pdfWeights)

process.h2ganalyzer.globalCounters.extend(['processedEvents'])

# PFIsolation photons
#process.load("HiggsAnalysis.HiggsTo2photons.pfIsolation_cff")

process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)

#################################################
#       Charged Hadron Subtraction Jets         #
#################################################

from RecoJets.JetProducers.ak5PFJets_cfi import *
process.ak5PFchsJets = ak5PFJets.clone()
process.ak5PFchsJets.src = 'pfNoPileUp'

#And the sequence below re-runs the pfNoPu on the fly
process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

# note pfPileUp modified according to JetMET's recommendations
process.pfPileUp.checkClosestZVertex = False
process.pfPileUp.Vertices = 'goodOfflinePrimaryVertices'

process.load("CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi")
process.pfNoPileUpSequence.insert(0, process.goodOfflinePrimaryVertices)

process.ak5PFchsL1Fastjet = cms.ESProducer(
    'L1FastjetCorrectionESProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK5PFchs'),
    srcRho      = cms.InputTag('kt6PFJets','rho')
)

process.ak5PFchsL2Relative = process.ak5CaloL2Relative.clone( algorithm = 'AK5PFchs' )
process.ak5PFchsL3Absolute     = process.ak5CaloL3Absolute.clone( algorithm = 'AK5PFchs' )
process.ak5PFchsResidual  = process.ak5CaloResidual.clone( algorithm = 'AK5PFchs' )

process.ak5PFchsL1FastL2L3 = cms.ESProducer(
  'JetCorrectionESChain',
  correctors =
  cms.vstring('ak5PFchsL1Fastjet','ak5PFchsL2Relative',
              'ak5PFchsL3Absolute')
)

process.ak5PFchsL1FastL2L3Residual = cms.ESProducer(
  'JetCorrectionESChain',
  correctors =
  cms.vstring('ak5PFchsL1Fastjet','ak5PFchsL2Relative',
              'ak5PFchsL3Absolute','ak5PFchsResidual')
)

#################################################
# Addition for Type1 pfMET corrections          #
#################################################

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

#################################################
# B-Tagging Modules                             #
#################################################

from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# for standard pfjets:
import RecoJets.JetAssociationProducers.ak5JTA_cff
process.newPFJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
process.newPFJetTracksAssociatorAtVertex.jets = "ak5PFJets"
process.newPFJetTracksAssociatorAtVertex.tracks = "generalTracks"

process.newPFImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
process.newPFImpactParameterTagInfos.jetTracks = "newPFJetTracksAssociatorAtVertex"
process.newPFTrackCountingHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
process.newPFTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )
#process.newPFTrackCountingHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
#process.newPFTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )
process.newPFJetProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
process.newPFJetProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )
#process.newPFJetBProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
#process.newPFJetBProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos") )

process.newPFSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
process.newPFSecondaryVertexTagInfos.trackIPTagInfos = "newPFImpactParameterTagInfos"
#process.newPFSimpleSecondaryVertexHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
#process.newPFSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFSecondaryVertexTagInfos") )
#process.newPFSimpleSecondaryVertexHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
#process.newPFSimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFSecondaryVertexTagInfos") )
process.newPFCombinedSecondaryVertexBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
process.newPFCombinedSecondaryVertexBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos"), cms.InputTag("newPFSecondaryVertexTagInfos") )
process.newPFCombinedSecondaryVertexMVABPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
process.newPFCombinedSecondaryVertexMVABPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFImpactParameterTagInfos"), cms.InputTag("newPFSecondaryVertexTagInfos") )

process.newPFJetTracksAssociator = cms.Sequence(
    process.newPFJetTracksAssociatorAtVertex
    )

process.newPFJetBtaggingIP = cms.Sequence(
    process.newPFImpactParameterTagInfos * (
       process.newPFTrackCountingHighEffBJetTags +
#       process.newPFTrackCountingHighPurBJetTags +
       process.newPFJetProbabilityBPFJetTags )
#       process.newPFJetBProbabilityBPFJetTags )
    )

process.newPFJetBtaggingSV = cms.Sequence(
    process.newPFImpactParameterTagInfos *
    process.newPFSecondaryVertexTagInfos * (
#       process.newPFSimpleSecondaryVertexHighEffBJetTags +
#       process.newPFSimpleSecondaryVertexHighPurBJetTags +
       process.newPFCombinedSecondaryVertexBPFJetTags +
       process.newPFCombinedSecondaryVertexMVABPFJetTags )
    )

process.newPFJetBtagging = cms.Sequence(
    process.newPFJetBtaggingIP +
    process.newPFJetBtaggingSV )

process.newPFBtaggingSequence = cms.Sequence(
    process.newPFJetTracksAssociator *
       process.newPFJetBtagging )



# for chs pfjets:
import RecoJets.JetAssociationProducers.ak5JTA_cff
process.newPFchsJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
process.newPFchsJetTracksAssociatorAtVertex.jets = "ak5PFchsJets"
process.newPFchsJetTracksAssociatorAtVertex.tracks = "generalTracks"

process.newPFchsImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
process.newPFchsImpactParameterTagInfos.jetTracks = "newPFchsJetTracksAssociatorAtVertex"
process.newPFchsTrackCountingHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
process.newPFchsTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsImpactParameterTagInfos") )
#process.newPFchsTrackCountingHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
#process.newPFchsTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsImpactParameterTagInfos") )
process.newPFchsJetProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
process.newPFchsJetProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsImpactParameterTagInfos") )
#process.newPFchsJetBProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
#process.newPFchsJetBProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsImpactParameterTagInfos") )

process.newPFchsSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
process.newPFchsSecondaryVertexTagInfos.trackIPTagInfos = "newPFchsImpactParameterTagInfos"
#process.newPFchsSimpleSecondaryVertexHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
#process.newPFchsSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsSecondaryVertexTagInfos") )
#process.newPFchsSimpleSecondaryVertexHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
#process.newPFchsSimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsSecondaryVertexTagInfos") )
process.newPFchsCombinedSecondaryVertexBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
process.newPFchsCombinedSecondaryVertexBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsImpactParameterTagInfos"), cms.InputTag("newPFchsSecondaryVertexTagInfos") )
process.newPFchsCombinedSecondaryVertexMVABPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
process.newPFchsCombinedSecondaryVertexMVABPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFchsImpactParameterTagInfos"), cms.InputTag("newPFchsSecondaryVertexTagInfos") )

process.newPFchsJetTracksAssociator = cms.Sequence(
    process.newPFchsJetTracksAssociatorAtVertex
    )

process.newPFchsJetBtaggingIP = cms.Sequence(
    process.newPFchsImpactParameterTagInfos * (
       process.newPFchsTrackCountingHighEffBJetTags +
#       process.newPFchsTrackCountingHighPurBJetTags +
       process.newPFchsJetProbabilityBPFJetTags )
#       process.newPFchsJetBProbabilityBPFJetTags )
    )

process.newPFchsJetBtaggingSV = cms.Sequence(
    process.newPFchsImpactParameterTagInfos *
    process.newPFchsSecondaryVertexTagInfos * (
#       process.newPFchsSimpleSecondaryVertexHighEffBJetTags +
#       process.newPFchsSimpleSecondaryVertexHighPurBJetTags +
       process.newPFchsCombinedSecondaryVertexBPFJetTags +
       process.newPFchsCombinedSecondaryVertexMVABPFJetTags )
    )

process.newPFchsJetBtagging = cms.Sequence(
    process.newPFchsJetBtaggingIP +
    process.newPFchsJetBtaggingSV )

process.newPFchsBtaggingSequence = cms.Sequence(
    process.newPFchsJetTracksAssociator *
       process.newPFchsJetBtagging )

#################################################
# Define path, first for AOD case then for RECO #
#################################################

#process.p11 = cms.Path(process.eventCounters*process.hcallLaserEvent2012Filter*process.ecalLaserCorrFilter*process.eventFilter1*process.pfNoPileUpSequence * process.pfParticleSelectionSequence * process.eleIsoSequence*process.ak5PFchsJets*process.producePFMETCorrections*process.newPFBtaggingSequence*process.newPFchsBtaggingSequence*process.eleRegressionEnergy * process.calibratedElectrons)
process.p11 = cms.Path(process.eventCounters*process.eventFilter1*process.pfNoPileUpSequence * process.pfParticleSelectionSequence * process.eleIsoSequence*process.ak5PFchsJets*process.producePFMETCorrections*process.newPFBtaggingSequence*process.newPFchsBtaggingSequence*process.eleRegressionEnergy * process.calibratedElectrons)
#process.p11 = cms.Path(process.eventCounters*process.eventFilter1* process.pfNoPileUpSequence * process.pfParticleSelectionSequence * process.eleIsoSequence*process.ak5PFchsJets*process.pfType1CorrectedMet  )

if (flagFastSim == 'OFF' or flagAOD == 'OFF'):
  process.p11 *= process.piZeroDiscriminators

#process.p11 *= (process.kt6PFJets* process.ak5PFJets* process.kt6PFJetsForRhoCorrection* process.h2ganalyzerPath)
process.p11 *= (process.h2ganalyzerPath)

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
  process.h2ganalyzer.doGenVertices = False
elif ((flagMC is 'ON' and flagFastSim is 'OFF') or flagPG is 'ON'):
  process.h2ganalyzer.doGenJet_algo1 = True
  process.h2ganalyzer.doGenJet_algo2 = True
  process.h2ganalyzer.doGenJet_algo3 = True
  process.h2ganalyzer.doGenParticles = True
  process.h2ganalyzer.doGenMet = True
  process.h2ganalyzer.doGenVertices = True
elif flagData is 'ON':
  process.h2ganalyzer.doPileup = False
  process.h2ganalyzer.doGenJet_algo1 = False
  process.h2ganalyzer.doGenJet_algo2 = False
  process.h2ganalyzer.doGenJet_algo3 = False
  process.h2ganalyzer.doGenParticles = False
  process.h2ganalyzer.doGenVertices = False
  process.h2ganalyzer.doGenMet = False

if ((flagMC is 'ON' or flagPG is 'ON') and flagAOD is 'OFF'):
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
  process.h2ganalyzer.EcalHitESColl = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES")
  process.h2ganalyzer.HcalHitsBEColl = cms.InputTag("hbhereco")
  process.h2ganalyzer.HcalHitsFColl = cms.InputTag("hfreco")
  process.h2ganalyzer.HcalHitsHoColl = cms.InputTag("horeco")
  process.h2ganalyzer.BarrelBasicClusterColl = cms.InputTag("")
  process.h2ganalyzer.BarrelBasicClusterShapeColl = cms.InputTag("multi5x5BasicClusters","multi5x5BarrelShapeAssoc")
  process.h2ganalyzer.JetTrackAssociationColl_algo3 = cms.InputTag("kt4JetTracksAssociatorAtVertex")

if (flagPG is 'OFF'):
  process.h2ganalyzer.doL1 = True
  process.h2ganalyzer.doHLT = True
else:
  process.h2ganalyzer.doL1 = False
  process.h2ganalyzer.doHLT = False
  process.h2ganalyzer.doJet_algoPF3 = False
  process.h2ganalyzer.doParticleGun = True

process.GlobalTag.globaltag = "START53_V9::All"
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("JetCorrectionsRecord"),
           tag = cms.string("JetCorrectorParametersCollection_Summer13_V1_MC_AK5PF"),
           connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
           label = cms.untracked.string('AK5PF')
          ),
  cms.PSet(record = cms.string("JetCorrectionsRecord"),
           tag = cms.string("JetCorrectorParametersCollection_Summer13_V1_MC_AK5PFchs"),
           connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
           label = cms.untracked.string('AK5PFchs')
          ),
  cms.PSet(record = cms.string("JetCorrectionsRecord"),
           tag = cms.string("JetCorrectorParametersCollection_Summer13_V1_MC_AK7PF"),
           connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
           label = cms.untracked.string('AK7PF')
          ),
)
process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag("TriggerResults","", hltLabel)
process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","", hltLabel)
