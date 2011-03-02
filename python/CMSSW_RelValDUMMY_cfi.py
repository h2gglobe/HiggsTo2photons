
import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
#'file:/local/path/to/local/file.root',
'/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0024/F2DF9FC5-E00F-E011-88AF-0030486790C0.root',
    ) );


