import FWCore.ParameterSet.Config as cms

#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
    '/store/relval/CMSSW_3_5_4/RelValH130GGgluonfusion/GEN-SIM-RECO/START3X_V24-v1/0004/D48A6E41-A82B-DF11-ABA9-0017313F02F2.root',
    '/store/relval/CMSSW_3_5_4/RelValH130GGgluonfusion/GEN-SIM-RECO/START3X_V24-v1/0004/C2EFAB03-2D2C-DF11-B914-002618943985.root',
    '/store/relval/CMSSW_3_5_4/RelValH130GGgluonfusion/GEN-SIM-RECO/START3X_V24-v1/0004/9E196A44-A62B-DF11-BBB5-00248C0BE012.root',
    '/store/relval/CMSSW_3_5_4/RelValH130GGgluonfusion/GEN-SIM-RECO/START3X_V24-v1/0004/92571321-A82B-DF11-BB02-001A92811742.root',
    '/store/relval/CMSSW_3_5_4/RelValH130GGgluonfusion/GEN-SIM-RECO/START3X_V24-v1/0004/309C04D2-A52B-DF11-83A7-001A92810AEC.root'
    ) );


