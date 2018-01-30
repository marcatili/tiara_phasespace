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

#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "HistoManager.hh"
#include "EventAction.hh"
//#include "G4RunManagerKernel.hh"
//#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"

#include "RunAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


TrackingAction::TrackingAction(RunAction* runaction):fRunAction(runaction)
//:G4UserTrackingAction()
{ }

//TrackingAction::~TrackingAction()
//{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

   	EventAction* evtAct ;
	G4int ParentID;

	// Get analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 

 	G4String procName = " ";

	if((aTrack->GetDefinition()->GetPDGEncoding()) == 22){

		if(aTrack->GetCreatorProcess()!=0) {procName =aTrack->GetCreatorProcess()->GetProcessName();}

		
		//G4cout<<"procName = "<<procName<<G4endl;   

		//G4cout<<" track id = "<< aTrack->GetTrackID()<<"  procName = "<<procName<<G4endl;   
		//fpTrackingManager->SetStoreTrajectory(true); //and then? see example extended/runAndEvent/RE01

		if (evtAct->IsInelastic(procName)) {
			//these are gamma prompt generated
			fRunAction->AddGammaPrompt(); //G4cout<<"gamma prompt found"<<G4endl;
			
			/*
  			const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
			evnb=evt->GetEventID();

			trknb=aTrack->GetTrackID();
			ParentID=aTrack->GetParentID();

			PGvertex=aTrack->GetVertexPosition();
			PGvertexDir=aTrack->GetVertexMomentumDirection();
			*/
   			PGEkinVtx=aTrack->GetVertexKineticEnergy();
		
			/*					
			PGTime	= aTrack->GetGlobalTime(); //time since the current event began	
			
			G4double PGx = PGvertex.x(), PGy = PGvertex.y(), PGz = PGvertex.z();
			G4double PGdx = PGvertexDir.x(), PGdy = PGvertexDir.y(), PGdz = PGvertexDir.z();


			//FILL TREE	
			analysisManager->FillNtupleIColumn(2, 0, evnb);
			analysisManager->FillNtupleIColumn(2, 1, trknb);

			analysisManager->FillNtupleDColumn(2, 2, PGx/cm);
			analysisManager->FillNtupleDColumn(2, 3, PGy/cm);
			analysisManager->FillNtupleDColumn(2, 4, PGz/cm);

			analysisManager->FillNtupleDColumn(2, 5, PGdx);
			analysisManager->FillNtupleDColumn(2, 6, PGdy);
			analysisManager->FillNtupleDColumn(2, 7, PGdz);

			analysisManager->FillNtupleDColumn(2, 8, PGEkinVtx/MeV);

			analysisManager->FillNtupleDColumn(2, 9, PGTime/ns);
			//analysisManager->FillNtupleDColumn(2, 6, dir);
			analysisManager->FillNtupleIColumn(2, 10, ParentID);
			*/

			//histos
			analysisManager->FillH1(1, PGEkinVtx/MeV);

			//G4cout<<PGEkinVtx/MeV<<G4endl;
			//analysisManager->AddNtupleRow(2); 

		}
		//i could store all gamma prompt created in the head with their ekin and vertex pos, but how should i connect them to the trigger?
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{;}





