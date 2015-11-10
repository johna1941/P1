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
// $Id: P1DetectorConstruction.cc 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file P1DetectorConstruction.cc
/// \brief Implementation of the P1DetectorConstruction class

#include "P1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4TriangularFacet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "P1SensitiveDetector.hh"
#include "G4SDManager.hh"

#include <fstream>

P1DetectorConstruction::P1DetectorConstruction()
: fFibreLV(0)
{ }

P1DetectorConstruction::~P1DetectorConstruction()
{ }

G4VPhysicalVolume* P1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // Materials
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* neoprene  = nist->FindOrBuildMaterial("G4_NEOPRENE"); // As an example, we'll be more specific closer to the time. 
  G4Material* lucite  = nist->FindOrBuildMaterial("G4_LUCITE");
  G4Material* water  = nist->FindOrBuildMaterial("G4_WATER");
  

  // World
  G4Box* solidWorld =
    new G4Box("World",        //its name
       1.*m, 1.*m, 1.*m);     //its size
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible()); //This means that when it sets the scale of the world it will ignore this. 
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  // Orb
  G4String name = "orb"; // Orb is simple - solid w/ radius. G4Sphere can be set as hollow w/ sectors/segments, but we've began simple. 
  G4VSolid* orb = new G4Orb(name,5.*cm);
  G4LogicalVolume* orb_lv = new G4LogicalVolume(orb,water,name); //(eg.) Neoprene, can be changed to something more suitable in the future. 
  new G4PVPlacement(0,G4ThreeVector(),orb_lv,name,logicWorld,0,false); // Orb one inside logical world

  // Scintillator
  name = "scintillator";
  G4VSolid* scint = new G4Orb(name,4.*cm); //Another orb, inside of the outer orb. r = 4cm cf. r = 5cm
//Geant4 is hierarchical, so placing one substance inside of another will displace the orginal. The mother displaces the daughter. This is more efficient than specifying a hollow sphere. 
  G4LogicalVolume* scint_lv = new G4LogicalVolume(scint,lucite,name);
 new G4PVPlacement(0,G4ThreeVector(),scint_lv,name,orb_lv,0,false); // Orb two inside of Orb one. 

// Fibre
name = "fibre";
G4VSolid* fibre = new G4Tubs(name,0.,0.05*cm,1.*um,0,360.*deg);
fFibreLV = new G4LogicalVolume(fibre,lucite,name);
new G4PVPlacement(0,G4ThreeVector(0.,0.,-3.9*cm),fFibreLV,name,scint_lv,0,false); // It's good practise to ask the code to check (when placing) that it doesn't overlap anything. To find out how to do this, look at the G4PVPlacement section; should be an additional argument.


  //always return the physical World
  return physWorld;
}

void P1DetectorConstruction::ConstructSDandField()
{
  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* fibreSD = new P1SensitiveDetector("Fibre");
  pSDman->AddNewDetector(fibreSD);
  fFibreLV->SetSensitiveDetector(fibreSD);
}

