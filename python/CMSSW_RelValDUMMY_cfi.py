
import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
'file:/hadoop/cms/phedex/store/data/Run2011A/Photon/AOD/05Jul2011ReReco-ECAL-v1/0000/0845B870-AFA7-E011-B24F-001A6478AB7C.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/00A587ED-A490-E011-A4CB-E0CB4E1A1183.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0A7C80B8-A590-E011-AD8B-E0CB4E19F9A2.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0AD88479-F790-E011-9BD9-E0CB4EA0A92E.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/20D9B651-A490-E011-A4F7-0030487CDAC2.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/2A8E7338-E290-E011-B9F6-E0CB4E29C4D9.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/3205E663-0691-E011-AA4C-001EC9D87221.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/3C7D9178-A690-E011-A9BF-E0CB4E29C4E9.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/407D77EA-A490-E011-8514-90E6BA0D09AA.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/480FAC7F-B490-E011-982E-E0CB4E19F9BD.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/4AF4444D-F390-E011-8130-485B39800B65.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/500849CE-A490-E011-8ABE-E0CB4E1A1191.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/70677204-A890-E011-B8FD-E0CB4E4408F7.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/723F010B-A790-E011-8AA7-E0CB4E19F9B9.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/78A399FF-A590-E011-ADFF-E0CB4E29C4F3.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/80130184-A690-E011-96AE-E0CB4E29C4C5.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/8253E475-F390-E011-AC2A-E0CB4E553644.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/8AB8442A-EF90-E011-B090-485B39800B65.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/8C621FFA-A790-E011-979C-485B39800BD2.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/8C7AC5E3-A590-E011-9C52-E0CB4EA0A936.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/909269C8-CE90-E011-A468-90E6BA0D09B2.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/92AA8EBF-A790-E011-ABF1-485B39800BC0.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/9E933311-A890-E011-BCF0-485B39800C17.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/A87A365C-A590-E011-92E3-E0CB4E1A114E.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/A8E24AF7-A690-E011-B3F0-E0CB4E19F96E.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/AA4A1801-A790-E011-A2BB-E0CB4EA0A8E0.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/BCFE0E33-A990-E011-B8C7-E0CB4E29C502.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/D6F95F0B-A790-E011-9D17-90E6BA442F15.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/D8391A28-DE90-E011-8CCD-485B39800BB0.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F4186041-A790-E011-8798-E0CB4E29C4C2.root',
'file:/hadoop/cms/phedex/store/mc/Summer11/VBF_HToGG_M-120_7TeV-powheg-pythia6/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FE8F05A6-A390-E011-85D9-00261834B5D2.root',
#'/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/161/217/1CABCC0B-E656-E011-9A1C-000423D996C8.root',

#'/store/mc/Spring11/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0012/84419D32-6550-E011-8980-00266CF2679C.root',
    #'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_1_2/RelValH130GGgluonfusion/GEN-SIM-RECO/START311_V2-v1/0019/00630515-FC44-E011-9757-001A92811748.root',
    #'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_11_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START311_V2-v1/0007/9E18EC03-294E-E011-BC4C-003048678E6E.root',
) );


