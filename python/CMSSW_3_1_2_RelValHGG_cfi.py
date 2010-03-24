import FWCore.ParameterSet.Config as cms

#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( ( 
        '/store/relval/CMSSW_3_1_2/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/FC6C2AF4-E278-DE11-B2D1-001D09F23A07.root',
        '/store/relval/CMSSW_3_1_2/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/CEA857A2-CC78-DE11-8073-000423D98EA8.root',
        '/store/relval/CMSSW_3_1_2/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/BA7034AD-CC78-DE11-96E0-001D09F251E0.root',
        '/store/relval/CMSSW_3_1_2/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/748489A8-CC78-DE11-991C-000423D99896.root',
        '/store/relval/CMSSW_3_1_2/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/104E25AC-CC78-DE11-AE55-001D09F2447F.root'
       ) );


