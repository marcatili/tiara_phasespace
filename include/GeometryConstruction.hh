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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef GeometryConstruction1_h
#define GeometryConstruction1_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4UserLimits.hh"
#include "GeometryConstructionMessenger.hh"

class G4VPhysicalVolume;
class GeometryConstructionMessenger;	

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GeometryConstruction : public G4VUserDetectorConstruction
{
  public:
    GeometryConstruction ();
   ~GeometryConstruction ();

    virtual void ConstructSDandField();
    virtual G4VPhysicalVolume* Construct();
    void SetMaxStep (G4double );
	
    G4double ShellThickness;
    void SetPhSPThickness(G4double val){ShellThickness= val;}; //set half thickness

  private:

	
    	GeometryConstructionMessenger* theMessenger;  

 	G4VPhysicalVolume* fworld_phys;

	G4LogicalVolume * head_log ;
	G4VPhysicalVolume* head_phys;

	G4LogicalVolume * skull_log ;
	G4VPhysicalVolume* skull_phys;

    //    G4VPhysicalVolume*  fDia_phys;
    //	G4LogicalVolume* Dia_log;

     //   G4VPhysicalVolume*  fAbs_phys;
     //	G4LogicalVolume* Abs_log;

        G4VPhysicalVolume*   PhSpShell_phys;
	G4LogicalVolume*  PhSpShell_log;


	G4UserLimits *fStepLimit;

};

    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

