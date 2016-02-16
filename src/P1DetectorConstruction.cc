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
#include "P1DetectorMessenger.hh"
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
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
// The above three are for the creation of ABS

#include <fstream>

P1DetectorConstruction::P1DetectorConstruction()
: fpDetectorMessenger(new P1DetectorMessenger(this))
//, fReflectivity(-1.) // (-1.) initialises it to -1, which is physically impossible. This is a good check to make sure that you've set it.
{ }

P1DetectorConstruction::~P1DetectorConstruction()
{
// delete fpDetectorMessenger; 
}

G4VPhysicalVolume* P1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;



  ///////////////////////////////////
  //////////// Materials ////////////
  //////////////////////////////////
  ////////// Construct ABS //////////
  G4String name, symbol;
  G4double density;
  G4int ncomponents, natoms;
  G4double fractionmass;
  G4UnitDefinition::BuildUnitsTable();

  // Carbon
  G4Element* C = nist->FindOrBuildElement("C");
  // Hydrogen
  G4Element* H = nist->FindOrBuildElement("H");
  // Nitrogen
  G4Element* N = nist->FindOrBuildElement("N");
  //Styrene
  density = 0.909*g/cm3;
  G4Material* styrene = new G4Material(name = "Styrene", density, ncomponents=2);
  styrene->AddElement(C, natoms=8); 
  styrene->AddElement(H, natoms=8);
  //1,3-Butadiene
  density = 0.6149*g/cm3; // At 25\degree (solid)
  G4Material* buta = new G4Material(name = "1,3-Butadiene", density, ncomponents=2);
  buta->AddElement(C, natoms=4);
  buta->AddElement(H, natoms=6);
  //Acrylonitrile
  density = 0.81*g/cm3;
  G4Material* acryl = new G4Material(name = "Acrylonitrile", density, ncomponents=3);
  acryl->AddElement(C, natoms=3); 
  acryl->AddElement(H, natoms=3);
  acryl->AddElement(N, natoms=1);

  // ABS
  density = 1.08*g/cm3; //1.06-1.08, according to wikipedia
  G4Material* ABS = new G4Material(name = "ABS", density, ncomponents=3); // Do not need to call it "G4_...", since this
  // is typically reserved for G4 library files. 
  ABS->AddMaterial(styrene, fractionmass=55*perCent); // 40-60%
  ABS->AddMaterial(buta, fractionmass=20*perCent); // 5-30%
  ABS->AddMaterial(acryl, fractionmass=25*perCent); //15-35%

  ////////// Construct LiqScint //////////
  // The atomic composition of our liquid scintillator is simply given as a
  // ratio of H:C (1.33), since there are no other molecules and its simply a 
  // percentage mixture of two hydrocarbons
  density = 1.136*g/cm3;
  G4Material* LS = new G4Material(name = "LS", density, ncomponents=2);
  LS->AddElement(C, natoms=3);
  LS->AddElement(H, natoms=4);





  // Materials
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  // These need to be subbed out for the true values.
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
// Create material properties table and add properties
  G4MaterialPropertiesTable* scint_mpt = new G4MaterialPropertiesTable();
  // Add to material properties table
  scint_mpt->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
  ->SetSpline(true);
  scint_mpt->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
  ->SetSpline(true);
  scint_mpt->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
  ->SetSpline(true);
  scint_mpt->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
  ->SetSpline(true);
  scint_mpt->AddConstProperty("SCINTILLATIONYIELD",0.1/MeV);
  scint_mpt->AddConstProperty("RESOLUTIONSCALE",1.0);
  scint_mpt->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  scint_mpt->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  scint_mpt->AddConstProperty("YIELDRATIO",0.8);
  G4cout << "Scint G4MaterialPropertiesTable\n"; scint_mpt->DumpTable();
  // Associate material properties table with the liquid scintillator material
  LS->SetMaterialPropertiesTable(scint_mpt);

  // Optical properties of the surface of the scintillator
  G4OpticalSurface* scint_surface = new G4OpticalSurface("scint-surface");
  scint_surface->SetType(dielectric_dielectric); // If both surfaces have refractive properties added, this will actually calculate reflection for us
  scint_surface->SetFinish(groundfrontpainted);
  scint_surface->SetModel(unified);
  G4cout << "scint_surface\n"; scint_surface->DumpInfo();
  // Create material properties table and add properties
  if (fReflectivity < 0.) {
    G4cout << "Reflectivity not set!" << G4endl;
    abort();
  }
  G4double reflectivity[nEntries]; for (auto& r: reflectivity) r = fReflectivity;
  G4MaterialPropertiesTable* mptForSkin = new G4MaterialPropertiesTable();  
  mptForSkin->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, nEntries)
  ->SetSpline(true);
  G4cout << "Skin G4MaterialPropertiesTable\n"; mptForSkin->DumpTable();
  // Associates the material properties with the surface of the liquid scintillator. 
  scint_surface->SetMaterialPropertiesTable(mptForSkin); 


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
  // G4String name = "orb"; // Orb is simple - solid w/ radius. G4Sphere can be set as hollow w/ sectors/segments, but we've began simple. 
  G4VSolid* orb = new G4Orb(name="orb",5.*cm);
  G4LogicalVolume* orb_lv = new G4LogicalVolume(orb,ABS,name);
  new G4PVPlacement(0,G4ThreeVector(),orb_lv,name,logicWorld,0,false); // Orb one inside logical world

  // Scintillator
  // name = "scintillator";
  G4VSolid* scint = new G4Orb(name="scintillator",4.*cm); //Another orb, inside of the outer orb. r = 4cm cf. r = 5cm
                                                          //Geant4 is hierarchical, so placing one substance inside of another will displace the orginal. The mother displaces the daughter. This is more efficient than specifying a hollow sphere.
  G4LogicalVolume* scint_lv = new G4LogicalVolume(scint,LS,name);
  new G4PVPlacement(0,G4ThreeVector(),scint_lv,name,orb_lv,0,false); // Orb two inside of Orb one.
                                                                     // Associate the optical surface
 
  ///////////////////                                                                  //  new G4LogicalSkinSurface("scint-surface", scint_lv, scint_surface);
  /// Central Orb ///
  ///////////////////
  
  // G4String name = "CentOrb"; // Orb is simple - solid w/ radius. G4Sphere can be set as hollow w/ sectors/segments, but we've began simple. 
  G4VSolid* CentOrb = new G4Orb(name="CentOrb",2.*cm);
  G4LogicalVolume* CentOrb_lv = new G4LogicalVolume(CentOrb,ABS,name); 
  new G4PVPlacement(0,G4ThreeVector(),CentOrb_lv,name,scint_lv,0,false);

  // Central Scintillator
  // name = "CentScintillator";
  G4VSolid* CentScint = new G4Orb(name="CentScintillator",1.9*cm); //Another orb, inside of the outer orb. r = 1.9 cm cf. r = 2 cm
                                                          //Geant4 is hierarchical, so placing one substance inside of another will displace the orginal. The mother displaces the daughter. This is more efficient than specifying a hollow sphere.
  G4LogicalVolume* CentScint_lv = new G4LogicalVolume(CentScint,LS,name);
  new G4PVPlacement(0,G4ThreeVector(),CentScint_lv,name,CentOrb_lv,0,false); // Orb two inside of Orb one.
                                                                     // Associate the optical surface
                                                                     //  new G4LogicalSkinSurface("scint-surface", scint_lv, scint_surface);




  // NOTE TO SC: To add an optical properties table, consult the code above, and the Geant4 application developers guide
  // Fibre1
  name = "fibre";
  G4VSolid* fibre = new G4Tubs(name,0.,5.*mm,1.*um,0,360.*deg);
  fFibreLV = new G4LogicalVolume(fibre,ABS,name);
  G4Transform3D transform = G4Translate3D(-2.25*cm,2.25*cm,-2.25*cm) * G4RotateZ3D(315.*deg) * G4RotateY3D(-135.*deg);
  fFibrePV = new G4PVPlacement(transform,fFibreLV,name,scint_lv,0,false,true);
  fFibre_axis = G4ThreeVector(0,0,1);

  // Fibre2
  name = "fibre2";
  G4VSolid* fibre2 = new G4Tubs(name,0.,0.05*cm,1.*um,0,360.*deg);
  fFibre2LV = new G4LogicalVolume(fibre2,ABS,name);
  G4Transform3D transform2 = G4Translate3D(2.3*cm,-2.3*cm,-2.3*cm) * G4RotateZ3D(315.*deg) * G4RotateY3D(135.*deg);
  fFibre2PV = new G4PVPlacement(transform2,fFibre2LV,name,scint_lv,0,false,true);
  fFibre2_axis = G4ThreeVector(0,1,0);

  // Fibre3
  name = "fibre3";
  G4VSolid* fibre3 = new G4Tubs(name,0.,0.05*cm,1.*um,0,360.*deg);
  fFibre3LV = new G4LogicalVolume(fibre3,ABS,name);
  G4Transform3D transform3 = G4Translate3D(-2.3*cm,-2.3*cm,2.3*cm) * G4RotateZ3D(45.*deg) * G4RotateY3D(-45.*deg);
  fFibre3PV = new G4PVPlacement(transform3,fFibre3LV,name,scint_lv,0,false,true);
  fFibre3_axis = G4ThreeVector(0,1,0);

  // Fibre4
  name = "fibre4";
  G4VSolid* fibre4 = new G4Tubs(name,0.,0.05*cm,1.*um,0,360.*deg);
  fFibre4LV = new G4LogicalVolume(fibre4,ABS,name);
  G4Transform3D transform4 = G4Translate3D(2.3*cm,2.3*cm,2.3*cm) * G4RotateZ3D(45.*deg) * G4RotateY3D(45.*deg);
  fFibre4PV = new G4PVPlacement(transform4,fFibre4LV,name,scint_lv,0,false,true);
  fFibre4_axis = G4ThreeVector(0,1,0);
  

  //always return the physical World
  return physWorld;
}

void P1DetectorConstruction::ConstructSDandField()
{
  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* fibreSD = new P1SensitiveDetector("Fibre");
  pSDman->AddNewDetector(fibreSD);
  fFibreLV->SetSensitiveDetector(fibreSD);
  fFibre2LV->SetSensitiveDetector(fibreSD);
  fFibre3LV->SetSensitiveDetector(fibreSD);
  fFibre4LV->SetSensitiveDetector(fibreSD);
}

