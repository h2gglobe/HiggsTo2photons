
import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
#'file:/local/path/to/local/file.root',
#'/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/161/217/1CABCC0B-E656-E011-9A1C-000423D996C8.root',
'/store/mc/Spring11/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0012/84419D32-6550-E011-8980-00266CF2679C.root',
    #'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_1_2/RelValH130GGgluonfusion/GEN-SIM-RECO/START311_V2-v1/0019/00630515-FC44-E011-9757-001A92811748.root',
    #'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_11_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START311_V2-v1/0007/9E18EC03-294E-E011-BC4C-003048678E6E.root',
) );


