import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# Run below command to print decay tree
# cmsRun python/cmssw_treedraw.py inputFiles=file:miniAOD.root

process = cms.Process('treedraw')

options = VarParsing.VarParsing('analysis')
options.parseArguments()

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),                                                                 
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(True),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(True),
                                   #status = cms.untracked.vint32( {1, 2, 4, 21, 22, 44, 52, 62, 23, 51, 52, 71} )
                                   )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles
    #    'file:../trash/miniAOD.root'
    )
)

process.p = cms.Path(process.printTree)
