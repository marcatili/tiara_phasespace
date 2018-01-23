#ifndef GeometryConstructionMessenger_h
#define GeometryConstructionMessenger_h 1
 
#include "G4UImessenger.hh"
#include "globals.hh"
 
class G4UIcommand;
class G4UIdirectory;
//class G4UIcmdWithoutParameter;
//class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
//class G4UIcmdWith3Vector;
//class G4UIcmdWith3VectorAndUnit;
//class G4UIcmdWithAnInteger;
//class G4UIcmdWithADouble;
//class G4UIcmdWithABool;
//class G4UIcmdWithoutParameter;
 
class GeometryConstruction;


 class GeometryConstructionMessenger: public G4UImessenger
 {
 public:
   GeometryConstructionMessenger(GeometryConstruction*);
   ~GeometryConstructionMessenger();
 

void SetNewValue(G4UIcommand *command, G4String newValues);
   //    Identifies the command which has been invoked by the user, extracts the
  //    parameters associated with that command (held in newValues), and uses
 //    these values with the appropriate member function of G4GeneralParticleSource.
//G4String GetCurrentValue(G4UIcommand *command);

private:
  GeometryConstruction *DetectorConstruction;

 private: //commands
   G4UIdirectory*              phspDir;
   G4UIdirectory*              detectorDirectory;
   G4UIcmdWithADoubleAndUnit*  thicknessCmd;	
};

#endif
