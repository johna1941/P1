#include "P1DetectorMessenger.hh"

#include "P1DetectorConstruction.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include <sstream>

P1DetectorMessenger::P1DetectorMessenger(P1DetectorConstruction * myDet)
:myDetector(myDet)
{
  fpP1CommandDirectory = new G4UIcommand("/PDMPhys/",this);
  fpP1CommandDirectory->SetGuidance("P1 detector control.");

  fpP1SetDirectory = new G4UIcommand("/PDMPhys/set/",this);
  fpP1SetDirectory->SetGuidance("Set commands.");

  fpReflectivityCommand = new G4UIcmdWithADouble("/PDMPhys/set/reflectivity",this);
  fpReflectivityCommand->SetGuidance("Define reflectivity of chmaber walls.");
}

P1DetectorMessenger::~P1DetectorMessenger () {
  delete fpReflectivityCommand;
  delete fpP1SetDirectory;
  delete fpP1CommandDirectory;
}

void P1DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if (command == fpReflectivityCommand)
  {
    myDetector->fReflectivity = fpReflectivityCommand->GetNewDoubleValue(newValues);
  }
}

