
import FWCore.ParameterSet.Config as cms


# file organized like this:
# COLLECTIONS
# CUTS
# DEBUG
# SELECTORS
# BOOLS
# OTHER

h2ganalyzer = cms.EDAnalyzer(
    "GlobeAnalyzer",
    RootFileName = cms.string('prova.root'),
    JobMaker = cms.string('jobmaker unknown'),
    globalCounters = cms.vstring(),
    
    # COLLECTIONS
    GeneratorColl = cms.InputTag("generator"),
    GenParticlesColl = cms.InputTag("genParticles"),
    GlobeReducedGendRMin = cms.double(0.3),
    
    pileupInfoCollection = cms.InputTag("addPileupInfo"),
    
    GenJetColl_algo1 = cms.InputTag("ak5GenJets"),
    GenJetColl_algo2 = cms.InputTag("ak7GenJets"),
    GenJetColl_algo3 = cms.InputTag("kt4GenJets"),
    
    SimHitList = cms.VInputTag(cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
                               cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),
                               cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
                               cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"),
                               cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"), 
                               cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"),
                               cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"),
                               cms.InputTag("g4SimHits","TrackerHitsTIDHighTof"),
                               cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"),
                               cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"), 
                               cms.InputTag("g4SimHits","TrackerHitsTECLowTof"),
                               cms.InputTag("g4SimHits","TrackerHitsTECHighTof")),
    SimHitPixBarrelLowColl = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
    SimHitPixEndcapLowColl = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
    SimHitTIBLowColl = cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"),
    SimHitTOBLowColl = cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"),
    SimHitTIDLowColl = cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"),
    SimHitTECLowColl = cms.InputTag("g4SimHits","TrackerHitsTECLowTof"),
    
    SimTrackColl = cms.InputTag("g4SimHits"),
    SimTrackPlusSimVertex = cms.bool(False),                   
    SimVertexColl = cms.InputTag("g4SimHits"),
    
    TPColl = cms.InputTag("mergedtruth","MergedTrackTruth"),
    TVColl = cms.InputTag("mergedtruth","MergedTrackTruth"),
    
    PixelRecHitsColl = cms.InputTag("siPixelRecHits"),
    RphiRecHitsColl = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),                       
    MatchedRecHitsColl = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    StereoRecHitsColl = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    
    EcalHitEBColl = cms.InputTag("reducedEcalRecHitsEB"),
    EcalHitEEColl = cms.InputTag("reducedEcalRecHitsEE"),
    EcalHitESColl = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
        
    HcalHitsBEColl = cms.InputTag("particleFlowRecHitHCAL", "Cleaned", "RECO"),
    HcalHitsFColl = cms.InputTag("particleFlowRecHitHCAL", "Cleaned", "RECO"), 
    HcalHitsHoColl = cms.InputTag("particleFlowRecHitHCAL", "Cleaned", "RECO"),
    
    CaloTowerColl = cms.InputTag("towerMaker"),
    
    BarrelBasicClusterColl = cms.InputTag("",""),
    EndcapBasicClusterColl = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),    
    BarrelBasicClusterShapeColl = cms.InputTag("",""),
    BarrelHybridClusterShapeColl = cms.InputTag("hybridSuperClusters","hybridShapeAssoc"),
    EndcapBasicClusterShapeColl = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapShapeAssoc"),
    BarrelSuperClusterColl = cms.InputTag("multi5x5SuperClusters","multi5x5BarrelSuperClusters"),
    BarrelHybridClusterColl = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    EndcapSuperClusterColl = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    HybridSuperClusterColl = cms.InputTag("correctedHybridSuperClusters"),
    
    TrackColl = cms.InputTag("generalTracks"),
    GsfTrackColl = cms.InputTag("electronGsfTracks"),
    # Temporary
    TrackColl2 = cms.InputTag("generalTracks"),
    TrackColl3 = cms.InputTag("electronGsfTracks"),
    AssocLabel = cms.string('TrackAssociatorByHits'),
    
    ElectronColl_std = cms.InputTag("gsfElectrons"),        
    eIDLabels        = cms.VInputTag(cms.InputTag("eidLoose"),
                                     cms.InputTag("eidTight")),
                                
    PhotonCollStd = cms.InputTag("photons"),
    # this will be from 420 ConvertedPhotonColl = cms.InputTag("allConversions"),
    ConvertedPhotonColl = cms.InputTag("trackerOnlyConversions"),
    MuonColl = cms.InputTag("muons"),
    
    JetColl_algo1 = cms.InputTag("ak5CaloJets"),
    JetColl_algo2 = cms.InputTag("ak7CaloJets"),
    JetColl_algo3 = cms.InputTag("kt4CaloJets"),
    JetColl_algoPF1 = cms.InputTag("ak5PFJets"),
    JetColl_algoPF2 = cms.InputTag("ak7PFJets"),
    JetColl_algoPF3 = cms.InputTag("kt4PFJets"),
    bcBColl = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    bcEColl = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),
    tkColl = cms.InputTag("generalTracks"),                                                 

    JetTrackAssociationColl_algo1 = cms.InputTag("ak5JetTracksAssociatorAtVertex"),    
    JetTrackAssociationColl_algo2 = cms.InputTag("ak7JetTracksAssociatorAtVertex"),
    JetTrackAssociationColl_algo3 = cms.InputTag(""),
    
    rhoCorrection = cms.InputTag("kt6PFJets","rho"),

    PFCandidateColl = cms.InputTag("particleFlow"),
    
    CaloMETColl = cms.InputTag("met"),
    TcMETColl = cms.InputTag("tcMet"),
    PFMETColl = cms.InputTag("pfMet"),

    BeamSpot = cms.InputTag("offlineBeamSpot"),
    VertexColl_std = cms.InputTag("offlinePrimaryVerticesWithBS"),
    VertexColl_nobs = cms.InputTag("offlinePrimaryVertices"),
    
    L1EtMiss = cms.InputTag("l1extraParticles","MET"),    
    L1CenJet = cms.InputTag("l1extraParticles","Central"),
    L1ForJet = cms.InputTag("l1extraParticles","Forward"),    
    L1TauJet = cms.InputTag("l1extraParticles","Tau"),
    L1EMNonIso = cms.InputTag("l1extraParticles","NonIsolated"),
    L1EMIso = cms.InputTag("l1extraParticles","Isolated"),
    L1Mu = cms.InputTag("l1extraParticles"),
    L1GtReadoutRecordTag = cms.InputTag("gtDigis"),
    #L1GtObjectMapRecordTag = cms.InputTag("hltL1GtObjectMap"),
    L1GtObjectMapRecordTag = cms.InputTag("L1GlobalTriggerReadoutRecord"),
    HLTParameters = cms.PSet(
    useSecondaryTrigger = cms.bool(False),
    PrimaryTriggerResultsTag   = cms.InputTag("TriggerResults"      ,"","HLT"),
    TriggerResultsTag          = cms.InputTag("hltTriggerSummaryAOD","","HLT")
    ),
    # Use this in case 2 different HLTs are stored, like 8e29 and 1e31
    # useSecondaryTrigger = cms.bool(True),
    # PrimaryTriggerResultsTag   = cms.InputTag("TriggerResults"      ,"","HLT8E29"),
    # SecondaryTriggerResultsTag = cms.InputTag("TriggerResults"      ,"","HLT"    ),
    # TriggerResultsTag          = cms.InputTag("hltTriggerSummaryAOD","","HLT8E29")),
    
    # CUTS
    GeneratorCuts = cms.PSet(EtCut = cms.double(10.0)),
    GenJetCuts = cms.PSet(EtCut = cms.double(-1.0)),
    
    SimHitCuts = cms.PSet(EnergyCut = cms.double(0.0)),
    SimTrackCuts = cms.PSet(EnergyCut = cms.double(0.0)),
    
    EcalHitCuts = cms.PSet(BarrelEnergyCut = cms.double(0.08),
                           EndcapEnergyCut = cms.double(0.24),
                           PreEnergyCut = cms.double(-999.0),
                           EcalMaxDR = cms.double(0.5),
                           KeepOutsideCone = cms.bool(True)),
    HcalHitsCuts = cms.PSet(HBHEEnergyCut = cms.double(0.35),
                            HFEnergyCut = cms.double(1.0),
                            HOEnergyCut = cms.double(0.7),
                            HcalMaxDR = cms.double(0.6),
                            KeepOutsideCone = cms.bool(True)),
    
    BasicClusterCuts = cms.PSet(EnergyCut = cms.double(-1.0)),
    SuperClusterCuts = cms.PSet(EnergyCut = cms.double(-1.0)),
    CaloTowerCuts = cms.PSet(EtCut = cms.double(1.0)),
    
    TrackCuts = cms.PSet(PtCut = cms.double(-1.0)),
    TPCuts = cms.PSet(tpLvpCut = cms.double(99999.0),
                      tpEtaCut = cms.double(2.5),
                      tpTvpCut = cms.double(99999.0),
                      tpPtCut = cms.double(1.0),
                      tpPdgidCut = cms.vint32(11, 22)),
    IsoCuts = cms.PSet(InnerCone = cms.double(0.015),
                       OuterCone = cms.double(0.3)),
    
    ElectronCuts = cms.PSet(EtCut = cms.double(0.0)),
    MuonCuts = cms.PSet(PtCut = cms.double(-1.0)),
    PhotonCuts = cms.PSet(EtCut = cms.double(0.0)),
    ConvertedPhotonCuts = cms.PSet(EtCut = cms.double(0.0)),
    JetCuts = cms.PSet(EnergyCut = cms.double(0.0)),

    # DEBUG
    Debug_Level = cms.int32(0),
    
    # SELECTORS
    SelectorNumber = cms.int32(0),
    SelectorCuts0 = cms.PSet(RecoPhotonNumber = cms.int32(1000),
                             RecoMuonNumber = cms.int32(1000),
                             RecoElectronNumber = cms.int32(1000),
                             MCLepPhotNumber = cms.int32(1000),
                             MCetcutel= cms.double(10.),
                             MCetcutmu= cms.double(10.),
                             MCetcutph= cms.double(10.),
                             MCetacutel= cms.double(2.7),
                             MCetacutmu= cms.double(2.7),
                             MCetacutph= cms.double(2.7),
                             nMCel = cms.int32(2000),
                             nMCmu = cms.int32(2000),
                             nMCph = cms.int32(2000),
                             nMCelmu = cms.int32(2),
                             nMCelph = cms.int32(2000),
                             nMCmuph = cms.int32(2000),
                             nMCelmuph = cms.int32(2000),
                             Recoetcutel= cms.double(10.),
                             Recoetcutmu= cms.double(10.),
                             Recoetcutph= cms.double(10.),
                             Recoetacutel= cms.double(2.7),
                             Recoetacutmu= cms.double(2.7),
                             Recoetacutph= cms.double(2.7),
                             nRecoel = cms.int32(2000),
                             nRecomu = cms.int32(2000),
                             nRecoph = cms.int32(2000),
                             nRecoelmu = cms.int32(2),
                             nRecoelph = cms.int32(2000),
                             nRecomuph = cms.int32(2000),
                             nRecoelmuph = cms.int32(2000)),
    
    # RUNNING MODULES
    doGenerator = cms.bool(False),
    doReducedGen = cms.bool(True),
    doGenParticles = cms.bool(True),
    doGenVertices = cms.bool(True),
    
    doGenJet_algo1 = cms.bool(True),
    doGenJet_algo2 = cms.bool(True),
    doGenJet_algo3 = cms.bool(True),
    
    doSimHits = cms.bool(False),
    doSimTracks = cms.bool(False),
    doSimTrackPlusSimVertex = cms.bool(False),
    
    doTrackingParticles = cms.bool(False),
    
    doTkRecHits = cms.bool(False),
    doEcalRecHits = cms.bool(True),
    
    doPreshowerHits = cms.bool(False),
    doEcal = cms.bool(True),
    doHcal = cms.bool(False),
    doHFHcal = cms.bool(False),
    
    doCaloTower = cms.bool(True),
    
    doL1 = cms.bool(True),
    doHLT = cms.bool(True),
    
    doTracks = cms.bool(True),
    doGsfTracks = cms.bool(True),

    doVertices_std = cms.bool(True),
    doVertices_nobs = cms.bool(True),
    
    doElectron_std = cms.bool(True),
    doMuon = cms.bool(True),
    doPhoton = cms.bool(True),
    doAllConversions = cms.bool(True),
    doLeptons = cms.bool(True),
    
    doJet_algo1 = cms.bool(True),
    doJet_algo2 = cms.bool(True),
    doJet_algo3 = cms.bool(True),
    doJet_algoPF1 = cms.bool(True),
    doJet_algoPF2 = cms.bool(True),
    doJet_algoPF3 = cms.bool(True),

    doPFCandidates = cms.bool(True),
    PFIsoOuterCone = cms.double(0.4),
    
    doMet = cms.bool(True),
    dotcMet = cms.bool(True),
    doPFMet = cms.bool(True),
    doHt = cms.bool(True),
    
    doRho = cms.bool(True),
    doPileup = cms.bool(True),
    
    doFastSim = cms.bool(False),
    doAodSim  = cms.bool(True),
    
    storeGsfTracksOnlyIfElectrons = cms.bool(True),

    #PAT
    
    electronTag = cms.untracked.InputTag("cleanPatElectrons"),
    photonTag   = cms.untracked.InputTag("cleanPatPhotons"),
    muonTag     = cms.untracked.InputTag("cleanPatMuons"),
    jetTag   = cms.untracked.InputTag("cleanPatJets"),
    
    # MISCELLANEOUS
    TrackAssociatorParameters = cms.PSet(muonMaxDistanceSigmaX = cms.double(0.0),
                                         muonMaxDistanceSigmaY = cms.double(0.0),
                                         CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
                                         dRHcal = cms.double(1.0),
                                         # bool   accountForTrajectoryChangeMuon = true
                                         # matching requirements 
                                         dREcal = cms.double(1.0),
                                         CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
                                         # association types
                                         useEcal = cms.bool(True),
                                         # preselection requirements in theta-phi space
                                         # allowed range: 
                                         #   dTheta = +/- dR
                                         #   dPhi = +/- dR  
                                         # (track trajectory changes are taken into account for muons)
                                         dREcalPreselection = cms.double(1.0),
                                         HORecHitCollectionLabel = cms.InputTag("horeco"),
                                         dRMuon = cms.double(9999.0),
                                         crossedEnergyType = cms.string('SinglePointAlongTrajectory'),
                                         propagateAllDirections = cms.bool(True),
                                         muonMaxDistanceX = cms.double(5.0),
                                         muonMaxDistanceY = cms.double(5.0),
                                         useHO = cms.bool(True), ## RecoHits
                                         
                                         accountForTrajectoryChangeCalo = cms.bool(False),
                                         DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
                                         EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
                                         dRHcalPreselection = cms.double(1.0),
                                         useMuon = cms.bool(False), ## RecoHits
                                         
                                         useCalo = cms.bool(False), ## CaloTowers
                                         
                                         # input tags
                                         EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                                         dRMuonPreselection = cms.double(0.2),

    usePreshower = cms.bool(False),
    dRPreshowerPreselection = cms.double(0.2),
    trajectoryUncertaintyTolerance = cms.double(-1.0),
                                         
                                         truthMatch = cms.bool(False), ## debugging information
                                         
                                         HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
                                         useHcal = cms.bool(True) ## RecoHits
                                         )
    )

