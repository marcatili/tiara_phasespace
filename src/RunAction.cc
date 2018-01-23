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

#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "G4Run.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction(), fHistoManager(0), TotPhSpTrgsPrim(0),TotPhSpTrgsSec(0),  TotEvents(0)
{
 fHistoManager = new HistoManager(); 

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fHistoManager;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }  

/*
  GCnt=0;
  PGCnt=0;
  NCnt=0;
*/
  TotPhSpTrgsPrim=0;
  TotPhSpTrgsSec=0;

//  DiaTrgsPrim=0;
//  DiaTrgsSec=0;
//  AbsTrgs=0;
  TotEvents=0;
  TotGammaPrompt=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{
  //save histograms
  //      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }  
	
  G4cout<< "___________________________________________________"<<G4endl; 
  G4cout<< "Number of generated gamma prompt gammas: "<< TotGammaPrompt <<G4endl; 
  G4cout<< "Number of prymary particles in phase space: "<< TotPhSpTrgsPrim <<G4endl; 
  G4cout<< "Number of secondary particles in phase space: "<< TotPhSpTrgsSec <<G4endl; 
 
  //G4cout<< "Number of trigger on absorber: "<< AbsTrgs<<G4endl; 
  //G4cout<< "Fraction of secondaries absorbed in converter: "<<G4float(AbsTrgs-DiaTrgsSec)/G4float(AbsTrgs)<<G4endl;
	
  G4cout<< "Total number of events: "<< TotEvents<<G4endl; 	
  
  G4cout<< "Percentage of primary gamma prompt in phase space: "<< G4float(TotGammaPrompt)*100./G4float(TotEvents) <<"%"<<G4endl; 
  //G4cout<< "Converter efficiency: "<< G4float(AbsTrgs)*100./G4float(TotEvents)<<"%"<<G4endl; 	
  G4cout<< "___________________________________________________"<<G4endl; 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
