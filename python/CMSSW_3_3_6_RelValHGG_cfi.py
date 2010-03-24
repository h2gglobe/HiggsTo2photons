import FWCore.ParameterSet.Config as cms

#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
            '/store/relval/CMSSW_3_3_6/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP3X_V8H-v1/0009/BE853EA6-9EE4-DE11-8BC6-002618943856.root',
        '/store/relval/CMSSW_3_3_6/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP3X_V8H-v1/0008/DEF97A82-38E4-DE11-9F2B-002618943966.root',
        '/store/relval/CMSSW_3_3_6/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP3X_V8H-v1/0008/DC464220-38E4-DE11-9565-002618943884.root',
        '/store/relval/CMSSW_3_3_6/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP3X_V8H-v1/0008/B4938C9B-39E4-DE11-940B-003048678B7C.root',
        '/store/relval/CMSSW_3_3_6/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP3X_V8H-v1/0008/A675D0B1-38E4-DE11-A546-003048678D6C.root',
        '/store/relval/CMSSW_3_3_6/RelValH130GGgluonfusion/GEN-SIM-RECO/STARTUP3X_V8H-v1/0008/608DD74C-38E4-DE11-9D78-00261894397D.root'
    ) );


