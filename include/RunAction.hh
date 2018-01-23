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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class HistoManager;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
   RunAction();
  ~RunAction();

   virtual void BeginOfRunAction(const G4Run*);
   virtual void EndOfRunAction(const G4Run*);

  private:
   HistoManager* fHistoManager;

//   G4long GCnt;
//   G4long PGCnt;
//   G4long NCnt;
     //G4long DiaTrgsPrim;
     //G4long DiaTrgsSec;
     //G4long AbsTrgs;
     G4long TotPhSpTrgsPrim;
     G4long TotPhSpTrgsSec;
     G4long TotEvents; 	
     G4long TotGammaPrompt; 	

   public:
   inline void AddPhSpTriggerPrimary(){TotPhSpTrgsPrim++;}	
   inline void AddPhSpTriggerSecondary(){TotPhSpTrgsSec++;}	
 //  inline void AddAbsTrigger(){AbsTrgs++;}	
   inline void AddEvent(){TotEvents++;}	
   inline void AddGammaPrompt(){TotGammaPrompt++;}	

/*
   inline void AddPGamma(){PGCnt ++;}
   inline void AddGamma(){GCnt ++;}
   inline void AddNeutron(){NCnt ++;}
*/		

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif



