import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.ParticleSelectors.pfCandsForIsolation_cff  import *
from CommonTools.ParticleFlow.Isolation.pfPhotonIsolation_cff import *
from CommonTools.ParticleFlow.Isolation.pfPhotonIsolationFromDeposits_cff import *

pfSelectedPhotons = cms.EDFilter(
    "GenericPFCandidateSelector",
    src = cms.InputTag("particleFlow"),
    cut = cms.string("pdgId()==22 && mva_nothing_gamma>0")
)

isoValPhotonWithCharged03 = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithCharged"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    vetos = cms.vstring('0.02', 'Threshold(1.0)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

isoValPhotonWithCharged04 = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithCharged"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.02', 'Threshold(1.0)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

isoValPhotonWithNeutral03 = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithNeutral"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'), # 0.3333,
    vetos = cms.vstring('Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

isoValPhotonWithNeutral04 = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithNeutral"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'), # 0.3333,
    vetos = cms.vstring('Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

isoValPhotonWithPhotons03noveto = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithPhotons"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    vetos = cms.vstring(),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

isoValPhotonWithPhotons04noveto = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithPhotons"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring(),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

isoValPhotonWithPhotons03 = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithPhotons"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    vetos = cms.vstring('EcalBarrel:0.045', 
                        #'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
                        'EcalBarrel:AbsThresholdFromTransverse(0.08)',
                        'EcalEndcaps:AbsThreshold(0.100)',
                        'EcalEndcaps:0.070'), 
                        #'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)'), 
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

isoValPhotonWithPhotons04 = cms.EDProducer(
    "CandIsolatorFromDeposits",
    deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("isoDepPhotonWithPhotons"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('EcalBarrel:0.045', 
                        #'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
                        'EcalBarrel:AbsThresholdFromTransverse(0.08)',
                        'EcalEndcaps:AbsThreshold(0.100)',
                        'EcalEndcaps:0.070'), 
                        #'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)'), 
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
    )

pfPhotonIsoDepositsSequence = cms.Sequence(
    isoDepPhotonWithCharged   +
    isoDepPhotonWithNeutral   +
    isoDepPhotonWithPhotons   
)

pfPhotonIsolationFromDepositsSequence03 = cms.Sequence(
    isoValPhotonWithCharged03  +
    isoValPhotonWithNeutral03  +
    isoValPhotonWithPhotons03  +
    isoValPhotonWithPhotons03noveto
)

pfPhotonIsolationFromDepositsSequence04 = cms.Sequence(
    isoValPhotonWithCharged04  +
    isoValPhotonWithNeutral04  +
    isoValPhotonWithPhotons04  +
    isoValPhotonWithPhotons04noveto
)

pfPhotonIsolationSequence = cms.Sequence(
    pfPhotonIsoDepositsSequence +
    pfPhotonIsolationFromDepositsSequence03 +
    pfPhotonIsolationFromDepositsSequence04
    )

pfBasedPhotonIsoSequence = cms.Sequence(
    pfCandsForIsolationSequence +
    pfSelectedPhotons +
    pfPhotonIsolationSequence
    ) 

