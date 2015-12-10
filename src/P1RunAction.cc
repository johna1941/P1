//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: P1RunAction.cc 89630 2015-04-23 12:11:28Z gcosmo $
//
/// \file P1RunAction.cc
/// \brief Implementation of the P1RunAction class

#include "P1RunAction.hh"
#include "P1PrimaryGeneratorAction.hh"
#include "P1DetectorConstruction.hh"
#include "P1Run.hh"

#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

P1RunAction::P1RunAction()
{}

P1RunAction::~P1RunAction()
{}

G4Run* P1RunAction::GenerateRun()
{
  return new P1Run; 
}

void P1RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void P1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  const P1Run* p1Run = static_cast<const P1Run*>(run);
  G4int numberOfPhotons = p1Run->GetPhotons();

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const P1PrimaryGeneratorAction* generatorAction
  = static_cast<const P1PrimaryGeneratorAction*>
  (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  //
  if (IsMaster()) {
    /* G4cout
    << "\n--------------------End of Global Run-----------------------";  
  G4cerr
  << "\n The run consists of " << nofEvents << " " << runCondition
  << "\n Number of photons reaching sensitive detector: "*/
  G4cerr << numberOfPhotons << G4endl;
  /*<< "\n------------------------------------------------------------"
  << G4endl;
  }
  else {
    G4cout
    << "\n--------------------End of Local Run------------------------";
*/  }


}
