
import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
#'/store/mc/Winter10/GluGluToHToZZTo2L2Nu_M-300_7TeV-powheg-pythia6/AODSIM/E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/0001/08F05D93-2118-E011-B3CA-00E08179185F.root',    
'/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0024/F2DF9FC5-E00F-E011-88AF-0030486790C0.root',
#'/store/mc/Winter10/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6/GEN-SIM-RECO/E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/0000/04A28052-B526-E011-9845-00151796D6E4.root',
#'/store/data/Run2010A/EG/RECO/v4/000/143/657/00110D8B-76AE-DF11-8BCF-003048D2C108.root',
#'/store/data/Run2010A/EG/RECO/v4/000/143/657/040DF055-8CAE-DF11-BE1C-001617E30CD4.root',
#'/store/data/Run2010A/EG/RECO/v4/000/143/657/04E50028-D6AE-DF11-9D1D-003048F024E0.root',
#'/store/data/Run2010A/EG/RECO/v4/000/143/657/04EBD522-8DAE-DF11-B5A4-003048F110BE.root',
#'/store/data/Run2010A/EG/RECO/v4/000/143/657/0652D10C-C1AE-DF11-ABEC-003048F11DE2.root',
#'/store/data/Run2010A/EG/RECO/v4/000/143/657/06C998C4-D9AE-DF11-9CB6-003048F118AC.root',
#'/store/data/Run2010A/EG/RECO/v4/000/143/657/080F8956-B0AE-DF11-93F1-0030487C8CB8.root',
    ) );


