// *******************************************************************
// * License and Disclaimer                                          *
// *                                                                 *
// * This software is copyright of Geant4 Associates International   *
// * Ltd (hereafter 'G4AI'). It is provided under the terms and      *
// * conditions described in the file 'LICENSE' included in the      *
// * software system.                                                *
// * Neither the authors of this software system nor G4AI make any   *
// * representation or warranty, express or implied, regarding this  *
// * software system or assume any liability for its use.            *
// * Please see the file 'LICENSE' for full disclaimer and the       *
// * limitation of liability.                                        *
// *******************************************************************
// $Id$
// John Allison  28th August 2015

#include "P1SensitiveDetector.hh"

#include "P1EventAction.hh"
#include "P1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"

P1SensitiveDetector::P1SensitiveDetector(const G4String& name)
: G4VSensitiveDetector(name)
{}

G4bool P1SensitiveDetector::ProcessHits(G4Step* step,
                                        G4TouchableHistory*)
{
  G4Track* track = step->GetTrack();
  const G4ParticleDefinition* pPDef = track->GetParticleDefinition();
  if (pPDef != G4OpticalPhoton::Definition()) {return true;}

  // It's an optical photon - kill it!!!
  track->SetTrackStatus(fStopAndKill);

  // It's an optical photon

  G4RunManager* runManager = G4RunManager::GetRunManager();
  const G4VUserDetectorConstruction* dc = runManager->GetUserDetectorConstruction();
  const P1DetectorConstruction* constdc = static_cast<const P1DetectorConstruction*>(dc);
  P1DetectorConstruction* p1dc = const_cast<P1DetectorConstruction*>(constdc);

  G4StepPoint* preSP = step->GetPreStepPoint();
  const G4TouchableHandle& preTH = preSP->GetTouchableHandle();
  G4VPhysicalVolume* prePV = preTH->GetVolume();
  G4ThreeVector axis;
  if (prePV == p1dc->fFibrePV) {
    axis = p1dc->fFibre_axis;
  } else {
    axis = p1dc->fFibre2_axis;
  }

  G4ThreeVector direction = track->GetMomentumDirection();
  if (direction * axis < 0.995) {
    // Too far from axis - don't count?
    return true;
  }

  const G4UserEventAction* ea = runManager->GetUserEventAction();
  const P1EventAction* constp1ea = static_cast<const P1EventAction*>(ea);
  P1EventAction* p1ea = const_cast<P1EventAction*>(constp1ea);

  p1ea->AddPhoton();

  //  G4double eDep = step->GetTotalEnergyDeposit();

  //  G4Track* track = step->GetTrack();
  //  const G4ParticleDefinition* pPDef = track->GetParticleDefinition();
  //  const G4String& partName = pPDef->GetParticleName();

  //  G4StepPoint* preSP = step->GetPreStepPoint();
  //  G4double ke = preSP->GetKineticEnergy();
  //  const G4TouchableHandle& preTH = preSP->GetTouchableHandle();
  //  G4VPhysicalVolume* prePV = preTH->GetVolume();
  //  G4int copyNo = prePV->GetCopyNo();

  //  G4StepPoint* postSP = step->GetPostStepPoint();
  //  const G4VProcess* postProcess = postSP->GetProcessDefinedStep();
  //  const G4String& postProcessname = postProcess->GetProcessName();

  //  G4cout << "Cell: " << copyNo << ", adding: " << eDep/keV << " keV" << G4endl;
  //  p1ea->AddEdep(copyNo,eDep);
  
  return true;
}
