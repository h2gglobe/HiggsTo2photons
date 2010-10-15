import FWCore.ParameterSet.Config as cms

#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
    '/store/relval/CMSSW_3_8_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START38_V9-v1/0021/BE2A2363-95BF-DF11-97C2-003048678B30.root',
    '/store/relval/CMSSW_3_8_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START38_V9-v1/0022/00B72CD4-F6BF-DF11-A1A2-003048678FA6.root',
    '/store/relval/CMSSW_3_8_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START38_V9-v1/0022/84CA81E4-EFBF-DF11-B3C2-00261894396E.root',
    '/store/relval/CMSSW_3_8_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START38_V9-v1/0022/94A499E8-EEBF-DF11-B94E-0018F3D09676.root',
    '/store/relval/CMSSW_3_8_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START38_V9-v1/0022/C23C5899-E1BF-DF11-B87A-001A92971B7E.root',
    '/store/relval/CMSSW_3_8_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START38_V9-v1/0022/C6D6F602-3BC0-DF11-A45B-003048678FB8.root'
    ) );


