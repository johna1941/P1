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
  ///////////////////////////////////
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
    G4double photonEnergyRI[] =
    { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
        2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
        2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
        2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
        2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
        3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
        3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
        3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV }; //32
    // Values from EJ301-em.xls data table. To be incorporated...
    G4double photonEnergySF[] =
    { 3.108*eV, 3.089*eV, 3.069*eV, 3.051*eV,
        3.039*eV, 3.032*eV, 3.022*eV, 3.014*eV,
        2.995*eV, 2.988*eV, 2.981*eV, 2.974*eV,
        2.967*eV, 2.963*eV, 2.960*eV, 2.956*eV,
        2.952*eV, 2.946*eV, 2.939*eV, 2.935*eV,
        2.932*eV, 2.928*eV, 2.925*eV, 2.922*eV,
        2.918*eV, 2.908*eV, 2.891*eV, 2.882*eV,
        2.872*eV, 2.868*eV, 2.865*eV, 2.858*eV,
        2.841*eV, 2.825*eV, 2.820*eV, 2.814*eV,
        2.803*eV, 2.794*eV, 2.781*eV, 2.771*eV,
        2.763*eV, 2.747*eV, 2.732*eV, 2.717*eV,
        2.702*eV, 2.673*eV, 2.645*eV, 2.617*eV,
        2.590*eV, 2.563*eV, 2.537*eV, 2.511*eV,
        2.486*eV, 2.462*eV, 2.438*eV, 2.414*eV,
        2.391*eV }; //57
    G4double scintilFast[] =
    { 0.052, 0.078, 0.119, 0.168,
        0.210, 0.255, 0.353, 0.455,
        0.727, 0.797, 0.853, 0.892,
        0.925, 0.936, 0.948, 0.958,
        0.967, 0.981, 0.991, 0.996,
        1.000, 0.997, 0.993, 0.986,
        0.969, 0.895, 0.799, 0.762,
        0.729, 0.724, 0.721, 0.711,
        0.708, 0.704, 0.701, 0.697,
        0.680, 0.662, 0.632, 0.598,
        0.554, 0.490, 0.438, 0.401,
        0.364, 0.312, 0.273, 0.238,
        0.208, 0.183, 0.158, 0.136,
        0.119, 0.104, 0.092, 0.08,
        0.07 }; //57
    G4double refractiveIndex[] =
    { 1.3435, 1.344,  1.3445, 1.345,  1.3455, //5
        1.346,  1.3465, 1.347,  1.3475, 1.348, //10
        1.3485, 1.3492, 1.35,   1.3505, 1.351, //15
        1.3518, 1.3522, 1.3530, 1.3535, 1.354, //20
        1.3545, 1.355,  1.3555, 1.356,  1.3568, //25
        1.3572, 1.358,  1.3585, 1.359,  1.3595, //30
        1.36,   1.3608}; // 32
    G4double absorption[] =  // Quoted to be 2.5-3m bulk absorption, we'll assume worst case.
    { 2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //5
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //10
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //15
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //20
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //25
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //30
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //35
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //40
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //45
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //50
        2.5*m, 2.5*m, 2.5*m, 2.5*m, 2.5*m, //55
        2.5*m, 2.5*m}; //57
    // Health check
    const G4int nEntriesSF = sizeof(photonEnergySF)/sizeof(G4double);
    const G4int nEntriesRI = sizeof(photonEnergyRI)/sizeof(G4double);
    assert(sizeof(refractiveIndex) == sizeof(photonEnergyRI));
    assert(sizeof(absorption) == sizeof(photonEnergySF));
    assert(sizeof(scintilFast) == sizeof(photonEnergySF));
    // Create material properties table and add properties
    G4MaterialPropertiesTable* scint_mpt = new G4MaterialPropertiesTable();
    // Add to material properties table
    scint_mpt->AddProperty("RINDEX",       photonEnergyRI, refractiveIndex, nEntriesRI)
    ->SetSpline(true);
    scint_mpt->AddProperty("ABSLENGTH",    photonEnergySF,   absorption,      nEntriesSF)
    ->SetSpline(true);
    scint_mpt->AddProperty("FASTCOMPONENT",photonEnergySF,   scintilFast,     nEntriesSF)
    ->SetSpline(true);
    scint_mpt->AddConstProperty("SCINTILLATIONYIELD",12000/MeV);
    scint_mpt->AddConstProperty("RESOLUTIONSCALE",1.0);
    scint_mpt->AddConstProperty("FASTTIMECONSTANT", 3.2*ns); // Given to be 3.2ns in EJ-301 PDF
    //  scint_mpt->AddConstProperty("SLOWTIMECONSTANT",32.3*ns); // "First three components; 3.2, 32.3, 270...?"
    //  scint_mpt->AddConstProperty("YIELDRATIO",0.8); // Relative strength of the fast vs. slow, i.e. 80% scintillations are fast.
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
    G4double photonEnergyRe[] =
    {   3.315*eV, 3.261*eV, 3.208*eV, 3.157*eV,
        3.108*eV, 3.060*eV, 3.014*eV, 2.925*eV,
        2.841*eV, 2.763*eV, 2.688*eV, 2.617*eV,
        2.550*eV, 2.486*eV, 2.426*eV, 2.368*eV,
        2.313*eV, 2.260*eV, 2.210*eV, 2.162*eV,
        2.116*eV, 2.072*eV };
    
    G4double reflectivity[] =
    {   0.7, 0.8, 0.87, 0.899,
        0.92, 0.934, 0.945, 0.955,
        0.9575, 0.96, 0.962, 0.9625,
        0.964, 0.964, 0.964, 0.965,
        0.965, 0.965, 0.965, 0.9645,
        0.963, 0.962 };
    
    const G4int nEntriesRe = sizeof(photonEnergyRe)/sizeof(G4double);
    assert(sizeof(reflectivity) == sizeof(photonEnergyRe));
    G4MaterialPropertiesTable* mptForSkin = new G4MaterialPropertiesTable();
    mptForSkin->AddProperty("REFLECTIVITY", photonEnergyRe, reflectivity, nEntriesRe) // Takes photonEnergy, nEntries from liquid scint above. Looks like this is instead of using rindex of a new surface
    ->SetSpline(true);
    G4cout << "Skin G4MaterialPropertiesTable\n"; mptForSkin->DumpTable();
    // Associates the material properties with the surface of the liquid scintillator. 
    scint_surface->SetMaterialPropertiesTable(mptForSkin);





  ///////////////////////////////////
  /////////// Build World ///////////
  ///////////////////////////////////

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
  G4VSolid* orb = new G4Orb(name="orb",2.5*cm);
  G4LogicalVolume* orb_lv = new G4LogicalVolume(orb,ABS,name);
  new G4PVPlacement(0,G4ThreeVector(),orb_lv,name,logicWorld,0,false); // Orb one inside logical world

  // Scintillator
  // name = "scintillator";
  G4VSolid* scint = new G4Orb(name="scintillator",2.*cm); //Another orb, inside of the outer orb. r = 4cm cf. r = 5cm
                                                          //Geant4 is hierarchical, so placing one substance inside of another will displace the orginal. The mother displaces the daughter. This is more efficient than specifying a hollow sphere.
  G4LogicalVolume* scint_lv = new G4LogicalVolume(scint,LS,name);
  new G4PVPlacement(0,G4ThreeVector(),scint_lv,name,orb_lv,0,false); // Orb two inside of Orb one.
  // Associate the optical surface
  new G4LogicalSkinSurface("scint-surface", scint_lv, scint_surface);



  // NOTE TO SC: To add an optical properties table, consult the code above, and the Geant4 application developers guide
  // Fibre1
  name = "fibre";
  G4VSolid* fibre = new G4Tubs(name,0.,2.5*mm,1.*um,0.,360.*deg);
  // G4Tubs(G4String name, G4double RMin, G4double RMax, G4double Dz, G4double SPhi, G4double DPhi)
  // RMin: inner radius, RMax: outer radius, Dz: half-length in z, SPhi: Starting phi in rad, DPhi: Angle of segment in rad
  fFibreLV = new G4LogicalVolume(fibre,ABS,name);
  G4Transform3D transform = G4Translate3D(0.,0.,1.982*cm);
  fFibrePV = new G4PVPlacement(transform,fFibreLV,name,scint_lv,0,false,true);
  fFibre_axis = G4ThreeVector(0,0,1);

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

