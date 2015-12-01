#ifndef P1DetectorMessenger_h
#define P1DetectorMessenger_h 1

//#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithADouble;

class P1DetectorConstruction;

class P1DetectorMessenger: public G4UImessenger
{
public:
  P1DetectorMessenger(P1DetectorConstruction * myDet);
  ~P1DetectorMessenger ();
  void SetNewValue(G4UIcommand * command,G4String newValues);
private:
  P1DetectorConstruction * myDetector;
  G4UIcommand* fpP1CommandDirectory;
  G4UIcommand* fpP1SetDirectory;
  G4UIcmdWithADouble* fpReflectivityCommand;
};

#endif

