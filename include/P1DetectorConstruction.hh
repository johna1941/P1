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
// $Id: P1DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file P1DetectorConstruction.hh
/// \brief Definition of the P1DetectorConstruction class

#ifndef P1DetectorConstruction_h
#define P1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4ThreeVector.hh"
#include "globals.hh"

class P1DetectorMessenger;
class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class P1DetectorConstruction : public G4VUserDetectorConstruction
{
  friend class P1SensitiveDetector;
  friend class P1DetectorMessenger;

public:

  P1DetectorConstruction();
  virtual ~P1DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  virtual void ConstructSDandField();

private:

  P1DetectorMessenger* fpDetectorMessenger;

  G4double fReflectivity;

  G4LogicalVolume* fFibreLV;
  G4LogicalVolume* fFibre2LV;
  G4LogicalVolume* fFibre3LV;
  G4LogicalVolume* fFibre4LV;
  G4VPhysicalVolume* fFibrePV;
  G4VPhysicalVolume* fFibre2PV;
  G4VPhysicalVolume* fFibre3PV;
  G4VPhysicalVolume* fFibre4PV;
  G4ThreeVector fFibre_axis;
  G4ThreeVector fFibre2_axis;
  G4ThreeVector fFibre3_axis;
  G4ThreeVector fFibre4_axis;
};

#endif
