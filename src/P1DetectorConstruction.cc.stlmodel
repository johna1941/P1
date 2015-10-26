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
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include <fstream>

P1DetectorConstruction::P1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

P1DetectorConstruction::~P1DetectorConstruction()
{ }

namespace  {
  std::ifstream stlfile("stlfile.stl");

  char c0;
  union {
    char ch[4];
    int n;
    float f;
  } u;
  char header[81];

  void get_header() {
    for (int i = 0; i < 80; ++i) {
      if (!stlfile.get(c0)) abort();
      header[i] = c0;
    }
    header[80] = '\0';
  }

  int get_int() {
    for (int i = 0; i < 4; ++i) {
      if (!stlfile.get(c0)) abort();
      u.ch[i] = c0;
    }
    return u.n;
  }

  float get_float() {
    for (int i = 0; i < 4; ++i) {
      if (!stlfile.get(c0)) abort();
      u.ch[i] = c0;
    }
    return u.f;
  }
}

G4VPhysicalVolume* P1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // Materials
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* neoprene  = nist->FindOrBuildMaterial("G4_NEOPRENE");

  // World
  G4Box* solidWorld =
    new G4Box("World",        //its name
       1.*m, 1.*m, 1.*m);     //its size
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  // Read stl file and construct G4TessellatedSolid.
  get_header();
  // Replace spaces by '-'.
  for (size_t i = 0; i < std::strlen(header); ++i) {
    if (header[i] == ' ') header[i] = '-';
  }
  G4cout << "Header: \"" << header <<"\"" << G4endl;
  int nFacets = get_int();
  G4cout << "No of facets: " << nFacets << G4endl;
  G4TessellatedSolid* ts = new G4TessellatedSolid(header);
  for (int iFacet = 0; iFacet < nFacets; ++iFacet) {
    G4Normal3D normal(get_float(),get_float(),get_float());
    G4Point3D vertex1(get_float()*mm,get_float()*mm,get_float()*mm);
    G4Point3D vertex2(get_float()*mm,get_float()*mm,get_float()*mm);
    G4Point3D vertex3(get_float()*mm,get_float()*mm,get_float()*mm);
    if (iFacet < 10) {  // Print first ten facets
      G4cout << "Normal: " << normal << G4endl;
      G4cout << "Vertex1: " << vertex1 << G4endl;
      G4cout << "Vertex2: " << vertex2 << G4endl;
      G4cout << "Vertex3: " << vertex3 << G4endl;
    }
    G4TriangularFacet* tf =
    new G4TriangularFacet(vertex1,vertex2,vertex3,ABSOLUTE);
    ts->AddFacet(tf);
    if (!stlfile.get(c0)) abort();  // Two dummy characters
    if (!stlfile.get(c0)) abort();
  }
  ts->SetSolidClosed(true);
  G4LogicalVolume* lv = new G4LogicalVolume(ts,neoprene,header);
  new G4PVPlacement(0,G4ThreeVector(),lv,header,logicWorld,0,false);

  // That should be it
  if (stlfile.get(c0)) abort();

  //always return the physical World
  return physWorld;
}
