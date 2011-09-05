import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi import *
from RecoEgamma.ElectronIdentification.electronIdSequence_cff import *

eIDFilter = cms.EDFilter("eIDFilter",
                         eIDCiCLabels = cms.VInputTag(cms.InputTag('eIDLooseForZSelection'),),
                         eIDVBTFLabels = cms.VInputTag(cms.InputTag('eIDVBTF95ForZSelection'),),
                         electronCollection = cms.InputTag("goodElectronsOver5"),
                         useOR = cms.bool(True),
                         nElectrons = cms.int32(1)
                         )

twoHEEPEleEIDFilter = eIDFilter.clone()
twoHEEPEleEIDFilter.eIDCiCLabels = cms.VInputTag(cms.InputTag('eIDTightForZSelection'),)
twoHEEPEleEIDFilter.eIDVBTFLabels = cms.VInputTag(cms.InputTag('eIDHEEPForestarSelection'),)
twoHEEPEleEIDFilter.eIDVBTFLabels = cms.VInputTag(cms.InputTag('eIDVBTF95ForZSelection'),)
twoHEEPEleEIDFilter.nElectrons = cms.int32(2)
twoHEEPEleEIDFilter.electronCollection = cms.InputTag("goodElectronsOver5")

twoEleEIDFilter = eIDFilter.clone()
twoEleEIDFilter.eIDCiCLabels = cms.VInputTag(cms.InputTag('eIDLooseForZSelection'),)
twoEleEIDFilter.eIDVBTFLabels = cms.VInputTag(cms.InputTag('eIDVBTF95ForZSelection'),)
twoEleEIDFilter.nElectrons = cms.int32(2)
twoEleEIDFilter.electronCollection = cms.InputTag("goodElectronsOver8")

eIDVBTF80ForZSelection = eidCutBasedExt.clone()
eIDVBTF80ForZSelection.electronIDType = 'robust'
eIDVBTF80ForZSelection.electronQuality = 'tight'
eIDVBTF80ForZSelection.src = cms.InputTag("goodElectronsOver5")

eIDTightForZSelection = eidTight.clone()
eIDTightForZSelection.src = cms.InputTag("goodElectronsOver5")

eIDVBTF95ForZSelection = eidCutBasedExt.clone()
eIDVBTF95ForZSelection.electronIDType = 'robust'
eIDVBTF95ForZSelection.electronQuality = 'loose'
eIDVBTF95ForZSelection.src = cms.InputTag("goodElectronsOver5")

eIDHEEPForestarSelection = eidCutBasedExt.clone()
eIDHEEPForestarSelection .electronIDType = 'robust'
eIDHEEPForestarSelection .electronQuality = 'highenergy'
eIDHEEPForestarSelection .src = cms.InputTag("goodElectronsOver5")

eIDLooseForZSelection = eidLoose.clone()
eIDLooseForZSelection.src = cms.InputTag("goodElectronsOver5")

electronIdentificationFilter = cms.Sequence(eIDVBTF95ForZSelection * eIDLooseForZSelection * eIDFilter)
twoElectronsIdentificationFilter = cms.Sequence(eIDVBTF95ForZSelection * eIDLooseForZSelection * twoEleEIDFilter)
twoHEEPElectronsIdentificationFilter = cms.Sequence(eIDVBTF80ForZSelection * eIDTightForZSelection *eIDHEEPForestarSelection * eIDVBTF95ForZSelection * twoHEEPEleEIDFilter)
