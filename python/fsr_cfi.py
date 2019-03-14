import FWCore.ParameterSet.Config as cms
from FSRPhotonRecovery.FSRPhotons.FSRphotonSequence_cff import addFSRphotonSequence
from PhysicsTools.NanoAOD.common_cff import Var, P4Vars

def addFSRphotonNano(process):
    addFSRphotonSequence(process, 'slimmedMuons', "FSRPhotonRecovery/FSRPhotons/data/PhotonMVAv9_BDTG800TreesDY.weights.xml")
    process.fsrPhotonTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("FSRRecovery", "selectedFSRphotons"),
        cut = cms.string(""), #we should not filter on cross linked collections
        name = cms.string("FSRPhoton"),
        doc  = cms.string("FSR photons"),
        singleton = cms.bool(False),
        extension = cms.bool(False),
        variables = cms.PSet(
            P4Vars,
            photonMVAIdValue = Var("userFloat('photonMVAIdValue')",float,doc="Photon MVA ID value"),
            FSRphotonMVAValue = Var("userFloat('FSRphotonMVAValue')",float,doc="Photon FSR MVA value"),
            PFphotonIso03 = Var("userFloat('PFphotonIso03')",float,doc="Isolation"),
        )
    )
    
    process.FSRphotonSequence += process.fsrPhotonTable
    process.fsr_step = cms.Path(process.FSRphotonSequence)
    process.schedule.insert(0, process.fsr_step)
