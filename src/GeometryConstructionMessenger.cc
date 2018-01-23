

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
//
//

//=====================================================================

#include "G4UIdirectory.hh"
//#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
//#include "G4UIcmdWith3Vector.hh"
//#include "G4UIcmdWith3VectorAndUnit.hh"
//#include "G4UIcmdWithAnInteger.hh"
//#include "G4UIcmdWithADouble.hh"
//#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"

#include "G4Tokenizer.hh"

#include "GeometryConstruction.hh"
#include "GeometryConstructionMessenger.hh"


GeometryConstructionMessenger::GeometryConstructionMessenger
(GeometryConstruction *detector):DetectorConstruction(detector)
{
 
   phspDir = new G4UIdirectory("/phasespace/");
   phspDir->SetGuidance("phase space control commands.");
 
   detectorDirectory = new G4UIdirectory("/phasespace/geometry/");
   detectorDirectory->SetGuidance("phase space geometry control sub-directory");
 
   thicknessCmd = new G4UIcmdWithADoubleAndUnit("/phasespace/geometry/Thickness",this);
   thicknessCmd->SetGuidance("set phase space thickness");
   thicknessCmd->SetParameterName("PhSpThickness",false,false);
   thicknessCmd->SetRange("PhSpThickness >= 0.");
   thicknessCmd->SetDefaultUnit("mm");
   thicknessCmd->SetUnitCategory("Length");
   thicknessCmd->AvailableForStates(G4State_PreInit ,G4State_Init, G4State_Idle);	
  // thicknessCmd->SetToBeBroadcasted(false);
 
     }

 GeometryConstructionMessenger::~GeometryConstructionMessenger()
 {
   delete thicknessCmd;
   delete detectorDirectory;
   delete phspDir;     
}


void GeometryConstructionMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{

  
if (command==thicknessCmd){

    DetectorConstruction->SetPhSPThickness(thicknessCmd->GetNewDoubleValue(newValues)); 

} else {
   G4cout << "No /phasespace/geometry command entered" << G4endl;
}

}




