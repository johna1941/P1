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
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
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
  G4Material* neoprene  = nist->FindOrBuildMaterial("G4_MUSCLE_SKELETAL_ICRP"); // As an example, we'll be more specific closer to the time. 
  G4Material* liq_scint  = nist->FindOrBuildMaterial("G4_AIR");  // Again, an example.
G4Material* PbBalloon = nist->FindOrBuildMaterial("G4_Pb");

  // For now give liq_scint some optical properties (from examples/extended/optical/OpNovice).
  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  G4double photonEnergy[] =
  { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
    2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
    2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
    2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
    2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
    3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
    3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
    3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
  G4double refractiveIndex1[] =
  { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
    1.346,  1.3465, 1.347,  1.3475, 1.348,
    1.3485, 1.3492, 1.35,   1.3505, 1.351,
    1.3518, 1.3522, 1.3530, 1.3535, 1.354,
    1.3545, 1.355,  1.3555, 1.356,  1.3568,
    1.3572, 1.358,  1.3585, 1.359,  1.3595,
    1.36,   1.3608};
  G4double absorption[] =
//  {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
//    15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
//    45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
//    52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
//    30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
//    17.500*m, 14.500*m };
  { 3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m,
    3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m,
    3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m,
    3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m, 3.*m,
    3.*m, 3.*m, 3.*m, 3.*m };
  G4double scintilFast[] =
  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00 };
  G4double scintilSlow[] =
  { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
    7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
    3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
    4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
    7.00, 6.00, 5.00, 4.00 };
  // Health check
  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));
  assert(sizeof(absorption) == sizeof(photonEnergy));
  assert(sizeof(scintilFast) == sizeof(photonEnergy));
  assert(sizeof(scintilSlow) == sizeof(photonEnergy));
  // Add to material properties table
  mpt->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
  ->SetSpline(true);
  mpt->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
  ->SetSpline(true);
  mpt->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
  ->SetSpline(true);
  mpt->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
  ->SetSpline(true);
  mpt->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  mpt->AddConstProperty("RESOLUTIONSCALE",1.0);
  mpt->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  mpt->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  mpt->AddConstProperty("YIELDRATIO",0.8);
  mpt->DumpTable();
  liq_scint->SetMaterialPropertiesTable(mpt);

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
  G4LogicalVolume* orb_lv = new G4LogicalVolume(orb,neoprene,name); //(eg.) Neoprene, can be changed to something more suitable in the future. 
  new G4PVPlacement(0,G4ThreeVector(),orb_lv,name,logicWorld,0,false); // Orb one inside logical world

  // Scintillator
  name = "scintillator";
  G4VSolid* scint = new G4Orb(name,4.2*cm); //Another orb, inside of the outer orb. r = 4cm cf. r = 5cm
//Geant4 is hierarchical, so placing one substance inside of another will displace the orginal. The mother displaces the daughter. This is more efficient than specifying a hollow sphere.
  G4LogicalVolume* scint_lv = new G4LogicalVolume(scint,liq_scint,name);
  new G4PVPlacement(0,G4ThreeVector(),scint_lv,name,orb_lv,0,false); // Orb two inside of Orb one.
  
/*
G4OpticalSurface* scint_surface = new G4OpticalSurface("scint-surface");
  scint_surface->SetType(dielectric_dielectric);
  scint_surface->SetFinish(polishedfrontpainted);
  scint_surface->SetModel(unified);
  new G4LogicalSkinSurface("scint-surface", scint_lv, scint_surface);
  scint_surface->DumpInfo();
*/

// Lead Balloon
name = "balloon";
G4VSolid* balloon = new G4Orb(name,4.1925*cm);
G4LogicalVolume* balloon_lv = new G4LogicalVolume(balloon,PbBalloon,name);
new G4PVPlacement(0,G4ThreeVector(),balloon_lv,name,scint_lv,0,false);

// Fibre1
name = "fibre";
G4VSolid* fibre = new G4Tubs(name,0.,0.12*cm,1.*nm,0,360.*deg);
fFibreLV = new G4LogicalVolume(fibre,liq_scint,name);
new G4PVPlacement(0,G4ThreeVector(0.,0.,-4.195*cm),fFibreLV,name,scint_lv,0,false); // It's good practise to ask the code to check (when placing) that it doesn't overlap anything. To find out how to do this, look at the G4PVPlacement section; should be an additional argument.

/*// Fibre2
name = "fibre2";
G4VSolid* fibre2 = new G4Tubs(name,0.,0.05*cm,1.*um,0,360.*deg);
fFibre2LV = new G4LogicalVolume(fibre2,liq_scint,name);
new G4PVPlacement(0,G4ThreeVector(0.,3.8*cm,0.),fFibre2LV,name,scint_lv,0,false);*/
// Need to make this rotated so that it's in the plane of the surface


  //always return the physical World
  return physWorld;
}

void P1DetectorConstruction::ConstructSDandField()
{
  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* fibreSD = new P1SensitiveDetector("Fibre");
  pSDman->AddNewDetector(fibreSD);
  fFibreLV->SetSensitiveDetector(fibreSD);
  // fFibre2LV->SetSensitiveDetector(fibreSD);
}

