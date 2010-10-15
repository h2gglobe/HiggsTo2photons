
import FWCore.ParameterSet.Config as cms


# file organized like this:
# COLLECTIONS
# CUTS
# DEBUG
# SELECTORS
# BOOLS
# OTHER

h2ganalyzer = cms.EDAnalyzer("GlobeAnalyzer",
    RootFileName = cms.string('prova.root'),
    
    # COLLECTIONS
    GeneratorColl = cms.InputTag("generator"),
    GenParticlesColl = cms.InputTag("genParticles"),
    GlobeReducedGendRMin = cms.double(0.3),
    
    GenJetColl_it5 = cms.InputTag("iterativeCone5GenJets"),
    GenJetColl_it7 = cms.InputTag("kt4GenJets"),
    GenJetColl_mid = cms.InputTag("ak5GenJets"),
    
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
    
    EcalHitEBColl = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    EcalHitEEColl = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    EcalHitESColl = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
        
    HcalHitsBEColl = cms.InputTag("hbhereco"),
    HcalHitsFColl = cms.InputTag("hfreco"), 
    HcalHitsHoColl = cms.InputTag("horeco"),
    
    CaloTowerColl = cms.InputTag("towerMaker"),
    
    BarrelBasicClusterColl = cms.InputTag("multi5x5BasicClusters","multi5x5BarrelBasicClusters"),
    EndcapBasicClusterColl = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),    
    BarrelBasicClusterShapeColl = cms.InputTag("multi5x5BasicClusters","multi5x5BarrelShapeAssoc"),
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
    
    PhotonCollStd = cms.InputTag("photons"),
    ConvertedPhotonColl = cms.InputTag("convertedPhotons"),
    
    MuonColl = cms.InputTag("muons"),
    
    JetColl_mid = cms.InputTag("ak5CaloJets"),
    JetColl_it5 = cms.InputTag("iterativeCone5CaloJets"),
    JetColl_it7 = cms.InputTag("kt4CaloJets"),
    JetColl_it5pf = cms.InputTag("iterativeCone5PFJets"),
    JetColl_sis5pf = cms.InputTag("ak5PFJets"),
    JetColl_kt4pf = cms.InputTag("kt4PFJets"),
    bcBColl = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    bcEColl = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),
    tkColl = cms.InputTag("generalTracks"),                                                 
    JetTrackAssociationColl_mid = cms.InputTag("ak5JetTracksAssociatorAtVertex"),    
    JetTrackAssociationColl_it5 = cms.InputTag("iterativeCone5JetTracksAssociatorAtVertex"),
    JetTrackAssociationColl_it7 = cms.InputTag("kt4JetTracksAssociatorAtVertex"),

    CaloMETColl = cms.InputTag("met"),
    TcMETColl = cms.InputTag("tcMet"),
    PFMETColl = cms.InputTag("pfMet"),
    
    VertexColl_std = cms.InputTag("offlinePrimaryVerticesWithBS"),
    VertexColl_pix = cms.InputTag("pixelVertices"),
    
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
    HLTParameters = cms.PSet(FullHLT = cms.bool(False),
#was                             PrimaryTriggerResultsTag   = cms.InputTag("TriggerResults"      ,"","HLT8E29"),
#was                             useSecondaryTrigger = cms.bool(True),
#was                             SecondaryTriggerResultsTag = cms.InputTag("TriggerResults"      ,"","HLT"    ),
#was                             TriggerResultsTag          = cms.InputTag("hltTriggerSummaryAOD","","HLT8E29")),
                             PrimaryTriggerResultsTag   = cms.InputTag("TriggerResults"      ,"","HLT"),
                             useSecondaryTrigger = cms.bool(False),
                             SecondaryTriggerResultsTag = cms.InputTag("TriggerResults"      ,"","HLT"    ),
                             TriggerResultsTag          = cms.InputTag("hltTriggerSummaryAOD","","HLT")),
                             #TriggerResultsTag = cms.InputTag("hltTriggerSummaryRAW","","HLT")),   # Uncomment this and change FullHLT to True if you want full old Trigger Info.
    
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
    JetCuts = cms.PSet(EnergyCut = cms.double(10.0)),
    
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
    
    doGenJet_mid = cms.bool(True),
    doGenJet_it5 = cms.bool(True),
    doGenJet_it7 = cms.bool(True),
    
    doSimHits = cms.bool(True),
    doSimTracks = cms.bool(True),
    doSimTrackPlusSimVertex = cms.bool(True),
    
    doTrackingParticles = cms.bool(True),

    doEcalRecHits = cms.bool(True),
    
    doPreshowerHits = cms.bool(True),
    doEcal = cms.bool(True),
    doHcal = cms.bool(True),
    doHFHcal = cms.bool(True),
    
    doCaloTower = cms.bool(True),
    
    doL1 = cms.bool(True),
    doHLT = cms.bool(True),
    
    doTracks = cms.bool(True),
    doGsfTracks = cms.bool(True),

    doVertices_std = cms.bool(True),
    doVertices_pix = cms.bool(True),
    doVtxCompat = cms.bool(True),
    
    doElectron_std = cms.bool(True),
    doMuon = cms.bool(True),
    doPhoton = cms.bool(True),
    doConvertedPhoton = cms.bool(False),
    doLeptons = cms.bool(True),
    
    doJet_mid = cms.bool(True),
    doJet_it5 = cms.bool(True),
    doJet_it7 = cms.bool(True),
    doJet_it5pf = cms.bool(True), 
    doJet_sis5pf = cms.bool(True),
    doJet_kt4pf = cms.bool(True), 
    
    doMet = cms.bool(True),
    dotcMet = cms.bool(True),
    doPFMet = cms.bool(True),
    doHt = cms.bool(True),
    
    doFastSim = cms.bool(False),
    doAodSim  = cms.bool(False),
    doEgammaSummer09Skim = cms.bool(False),

    storeGsfTracksOnlyIfElectrons = cms.bool(True),
    
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
                                         ),
    
    
## Name of the Triggers from CMSSW_3_1_2/src/HLTrigger/Configuration/python/HLT_1E31_cff.py
MuonHLTLabels = cms.VInputTag(
    cms.InputTag("hltL1MuOpenL1Filtered0","","HLT"), #              HLT_L1MuOpen
    cms.InputTag("hltL1MuL1Filtered0","","HLT"), #                  HLT_L1Mu
    cms.InputTag("hltL1Mu20HQL1Filtered20","","HLT"), #             HLT_L1Mu20HQ
    cms.InputTag("hltL1Mu30L1Filtered30","","HLT"), #               HLT_L1Mu30
    cms.InputTag("hltL2Mu11L2Filtered11","","HLT"), #               HLT_L2Mu11
    cms.InputTag("hltSingleMuIsoL3IsoFiltered9","","HLT"), #        HLT_IsoMu9
    cms.InputTag("hltSingleMu5L3Filtered5","","HLT"), #             HLT_Mu5
    cms.InputTag("hltSingleMu9L3Filtered9","","HLT"), #             HLT_Mu9
    cms.InputTag("hltSingleMu11L3Filtered11","","HLT"), #           HLT_Mu11
    cms.InputTag("hltSingleMu15L3PreFiltered15","","HLT"), #        HLT_Mu15
    cms.InputTag("hltDoubleMuLevel1PathL1OpenFiltered","","HLT"), # HLT_L1DoubleMuOpen
    cms.InputTag("hltDiMuonL3PreFiltered0","","HLT"), #             HLT_DoubleMu0
    cms.InputTag("hltDiMuonL3PreFiltered","","HLT")  #              HLT_DoubleMu3
    ),
ElectronHLTLabels = cms.VInputTag(
    cms.InputTag("hltPreL1SingleEG5","","HLT"), #                                                   HLT_L1SingleEG5
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt10PixelMatchFilter","","HLT"), #              HLT_Ele10_SW_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15PixelMatchFilter","","HLT"), #              HLT_Ele15_SW_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15EleIdDphiFilter ","","HLT"), #              HLT_Ele15_SW_EleId_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15LTITrackIsolFilter","","HLT"), #            HLT_Ele15_SW_LooseTrackIso_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronSiStripEt15PixelMatchFilter","","HLT"), #       HLT_Ele15_SiStrip_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15LTIESDoubleSC15","","HLT"), #               HLT_Ele15_SC15_SW_LooseTrackIso_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15EleIdESDoubleSC15","","HLT"), #             HLT_Ele15_SC15_SW_EleId_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt20PixelMatchFilter","","HLT"), #              HLT_Ele20_SW_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronSiStripEt20PixelMatchFilter","","HLT"), #       HLT_Ele20_SiStrip_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt20ESDoubleSC15","","HLT"), #                  HLT_Ele20_SC15_SW_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15EtFilterESet25","","HLT"), #                HLT_Ele25_SW_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15EleIdTrackIsolFilterESet25LTI","","HLT"), # HLT_Ele25_SW_EleId_LooseTrackIso_L1R
    cms.InputTag("hltL1NonIsoDoubleElectronEt5JpsiPMMassFilter","","HLT"), #                        HLT_DoubleEle5_SW_Jpsi_L1R
    cms.InputTag("hltL1NonIsoDoubleElectronEt5UpsPMMassFilter","","HLT"), #                         HLT_DoubleEle5_SW_Upsilon_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoDoubleElectronEt10PixelMatchFilter","","HLT")  #              HLT_DoubleEle10_SW_L1R
    ),
PhotonHLTLabels = cms.VInputTag(
    cms.InputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter","","HLT"), # HLT_Photon10_L1R
    cms.InputTag("hltL1NonIsoHLTLEITISinglePhotonEt10TrackIsolFilter","","HLT"), # HLT_Photon10_LooseEcalIso_TrackIso_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter","","HLT"), # HLT_Photon15_L1R
    cms.InputTag("hltL1NonIsoHLTLEITISinglePhotonEt20TrackIsolFilter","","HLT"), # HLT_Photon20_LooseEcalIso_TrackIso_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt25HcalIsolFilter","","HLT"), # HLT_Photon25_L1R
    cms.InputTag("hltL1NonIsoHLTLEITISinglePhotonEt25TrackIsolFilter","","HLT"), # HLT_Photon25_LooseEcalIso_TrackIso_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt30HcalIsolFilter","","HLT"), # HLT_Photon30_L1R_1E31
    cms.InputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt30EtFilterESet70","","HLT"), # HLT_Photon70_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoDoublePhotonEt10HcalIsolFilter","","HLT"), # HLT_DoublePhoton10_L1R
    cms.InputTag("hltL1NonIsoHLTNonIsoDoublePhotonEt15HcalIsolFilter","","HLT"), # HLT_DoublePhoton15_L1R
    cms.InputTag("hltL1NonIsoHLTVLEIDoublePhotonEt15HcalIsolFilter"  ,"","HLT")  # HLT_DoublePhoton15_VeryLooseEcalIso_L1R
    )

)
