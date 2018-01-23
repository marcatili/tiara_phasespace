#include "EventAction.hh"
#include "G4Event.hh"
#include "HodoscopeHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "HistoManager.hh"
#include "RunAction.hh"

#include "G4VTrajectory.hh"
#include "G4RunManager.hh"

EventAction::EventAction(RunAction* runaction):fRunAction(runaction)
//:fHistoManager(0)
{

  //  fHistoManager = new HistoManager(); 


  //diaHCID = -1;
 // absHCID = -1;
   phspHCID = -1;
}

EventAction::~EventAction()
{
 
 // delete fHistoManager;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  	G4SDManager * SDman = G4SDManager::GetSDMpointer();

  	if((phspHCID<0)){
  		G4String colNam;
   		phspHCID = SDman->GetCollectionID(colNam="phasespace/HodoscopehitsCollection");

		//G4cout << ">>> COLLECTION ID " << phspHCID << G4endl;		   	
	} 
	

}


void EventAction::EndOfEventAction(const G4Event* evt)
{
	//G4cout << ">>> Event " << evt->GetEventID() << G4endl;
	fRunAction->AddEvent();
	

	if((phspHCID<0)) { 

		G4cout << "one of HCID<0!!!!!!!!!" << G4endl; 
		return;
	} 

	G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
	HodoscopeHitsCollection* phspHC = 0;
	//HodoscopeHitsCollection* absHC = 0;


	if(HCE){
	   	phspHC = (HodoscopeHitsCollection*)(HCE->GetHC(phspHCID));
	   	//absHC = (HodoscopeHitsCollection*)(HCE->GetHC(absHCID));

	} else {
		G4cout<<"WARNING : no hit collection for event " << evt->GetEventID() << G4endl;
		return;
	}


	if ((!phspHC) ) {
        	G4cout<<"WARNING : Some of hits collections of this event not found." << G4endl; 
        	return;
	}   


	// Get analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 	//analysisManager->FillNtupleIColumn(2, 0, evt->GetEventID() );

	///////////////////////////////
	//phase space
	///////////////////////////////

	G4long EvtNb=0;
	G4double ekin=0.;
	G4double x;
	G4double y;
	G4double z;
	G4double dx;
	G4double dy;
	G4double dz;
	G4double time;
	G4int PID;
	G4int ParentID;
	G4int TrackID;
	G4int VertexVolumeTag;
	G4int ProcessNameTag;
	G4String TrigPartName;
	G4String ProcessName;
	G4String VertexVolumeName;
	G4int IsGammaPrompt;

	/*
	G4double PGx;
	G4double PGy;
	G4double PGz;
	G4double PGdx;
	G4double PGdy;
	G4double PGdz;
	*/
  	
	//G4ThreeVector TrackVertexPos;
	//G4ThreeVector TrackVertexMom;
        //G4double VertexEkin;
	//G4double PGtime;

   	if ( (phspHC->entries())>0 ) {

		

		//G4cout << ">>> HIT collection of event " << evt->GetEventID() << G4endl;
		//get deposited energy and trigger time

		int n_hit = phspHC->entries();
		EvtNb=evt->GetEventID();

		for ( int i = 0; i < n_hit; i++){//loop on collection hits
			

			//initialise variables
			VertexVolumeTag=0;
			ProcessNameTag=0;
			IsGammaPrompt=0;

			//all particles in phase space

			//still it should be checked that:
			//	- the particle is entering the volume for the first time!!
			//      - it is not a primary proton			

			time   = (*phspHC)[i]->GetTime();
			ekin = (*phspHC)[i]->GetEnergy();

			x=((*phspHC)[i]->GetPartPosition()).x();
			y=((*phspHC)[i]->GetPartPosition()).y();
			z=((*phspHC)[i]->GetPartPosition()).z();


			dx=((*phspHC)[i]->GetPartDirection()).x();
			dy=((*phspHC)[i]->GetPartDirection()).y();
			dz=((*phspHC)[i]->GetPartDirection()).z();

			PID=(*phspHC)[i]->GetParticle()->GetPDGEncoding();
			ParentID=(*phspHC)[i]->GetPID();

			VertexVolumeName=(*phspHC)[i]->GetVertexVolumeName();
			//tag for vertex volume name
			if (VertexVolumeName=="world_L") VertexVolumeTag=1;
			if (VertexVolumeName=="head_log") VertexVolumeTag=2;
			if (VertexVolumeName=="skull_log") VertexVolumeTag=3;

			/*
			//infos at vertex
			TrackVertexPos =(*phspHC)[i]->GetTrackVertex();
			vx=TrackVertexPos.x(); vy=TrackVertexPos.y(); vz=TrackVertexPos.z();
						
			TrackVertexMom=	(*phspHC)[i]->GetTrackVertexMomentum();
			vdx=TrackVertexMom.x(); vdy=TrackVertexMom.y(); vdz=TrackVertexMom.z();

		        VertexEkin= (*phspHC)[i]->GetKinEnergy();
			*/

			TrackID= (*phspHC)[0]->GetTID(); 
			//The primary is the track with ID 1 ( and with no parent)
			//primary (has ParentID == 0)
			//if (TrackID!=1 ) ProcessName=(*phspHC)[i]->GetProcess(); 
			//	else ProcessName="primary"; 

			if ((TrackID==1) && (ParentID==0)) ProcessName="primary";
				else ProcessName=(*phspHC)[i]->GetProcess();  

			//tag for process name generating triggering particle
			if (ProcessName=="primary") ProcessNameTag=1;
			if (ProcessName=="compt") ProcessNameTag=2;
			if (ProcessName=="conv") ProcessNameTag=3;
			if (ProcessName=="annihil") ProcessNameTag=4;
			if (ProcessName=="hIoni") ProcessNameTag=5;
			if (ProcessName=="eIoni") ProcessNameTag=6;	
			
			
			//gamma	
			//triggering particle is photon -> look if it is prompt or scattered
			if ((PID==22)) {
				if(IsInelastic(ProcessName)){
					IsGammaPrompt=1; 
					ProcessNameTag=7;
					/*
					TrackVertexPos =(*phspHC)[i]->GetTrackVertex();
					PGx=TrackVertexPos.x(); PGy=TrackVertexPos.y(); PGz=TrackVertexPos.z();
						
					TrackVertexMom=	(*phspHC)[i]->GetTrackVertexMomentum();
					PGdx=TrackVertexMom.x(); PGdy=TrackVertexMom.y(); PGdz=TrackVertexMom.z();

					PGEkin= (*phspHC)[i]->GetKinEnergy();
					PGtime=(*phspHC)[i]->
					*/
	
				}else {  //if scattered look for vertex of parent track 
					//(in principle should be done until a gamma prompt is found- i.e. multiple compton scattering )
					//use this example: http://geant4.slac.stanford.edu/Tips/event/3.html
					IsGammaPrompt=0;	
					
					}	
			}	

			//print variables	
			/*
			G4cout << "------phase space Event # " <<EvtNb << G4endl;
			G4cout << "Particle ID " << PID << G4endl;
			G4cout << "Track ID " << TrackID << G4endl;
			G4cout << "Parent ID " << ParentID << G4endl;
			G4cout << "VertexVolumeTag " << VertexVolumeTag << G4endl;
			G4cout << "ProcessNameTag " << ProcessNameTag << G4endl;
			G4cout << "IsGammaPrompt " << IsGammaPrompt << G4endl;
			G4cout << "hit time " << (time/s)*1.E12<< " [ps] " <<G4endl;
			G4cout << "ekin " << ekin/MeV << " [MeV] " <<G4endl;
			G4cout << "x " << x/cm << " [cm] " <<G4endl;
			G4cout << "y " << y/cm << " [cm] " << G4endl;
			G4cout << "z " << z/cm << " [cm] " << G4endl;
			G4cout << "dx " << dx/cm << G4endl;
			G4cout << "dy " << dy/cm << G4endl;
			G4cout << "dz " << dz/cm << G4endl;
			*/			
			G4double r2=y/cm*y/cm+z/cm*z/cm;
			if ((PID!=2212) || ((PID==2212)&&(r2>0.01))){//fill tree


				//stille there may be protons scored multiple times
				//(because proton tracks are not killed in ProcessHits)
								
				//fill the tree		
				analysisManager->FillNtupleIColumn(1, 0, EvtNb);
				analysisManager->FillNtupleDColumn(1, 1, time/ns );
				analysisManager->FillNtupleDColumn(1, 2, ekin/MeV);
				analysisManager->FillNtupleDColumn(1, 3, x/cm);
				analysisManager->FillNtupleDColumn(1, 4, y/cm);
				analysisManager->FillNtupleDColumn(1, 5, z/cm);
				analysisManager->FillNtupleDColumn(1, 6, dx);//check if there are units!!!!
				analysisManager->FillNtupleDColumn(1, 7, dy);
				analysisManager->FillNtupleDColumn(1, 8, dz);
				analysisManager->FillNtupleIColumn(1, 9, PID);
				analysisManager->FillNtupleIColumn(1, 10, TrackID);
				analysisManager->FillNtupleIColumn(1, 11, VertexVolumeTag);
				analysisManager->FillNtupleIColumn(1, 12, ProcessNameTag);
				analysisManager->FillNtupleIColumn(1, 13, IsGammaPrompt);
				analysisManager->FillNtupleIColumn(1, 14, ParentID);


				//histos
				analysisManager->FillH1(2, time/ns);

				analysisManager->AddNtupleRow(1); 

			}

			//this a phase space trigger from prymary	
			if ((TrackID==1) && (ParentID==0))fRunAction->AddPhSpTriggerPrimary();  
			else fRunAction->AddPhSpTriggerSecondary();  //this a phase space trigger from secondary


		}

	}

}

//ionInelastic alphaInelastic tInelastic dInelastic neutronInelastic protonInelastic He3Inelastic
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
G4bool EventAction::IsInelastic(G4String process)
{
	//ok for g4.9.5
	//if(process =="IonInelastic" || process =="AlphaInelastic" || process=="ProtonInelastic" || process=="NeutronInelastic" || process=="DeuteronInelastic" || process=="TritonInelastic" || process=="He3Inelastic") return true;

	// wrong!!
	//if((process =="ionInelastic") || (process =="alphaInelastic") || (process=="protonInelastic") || (process=="neutronInelastic") || (process=="deuteronInelastic") || (process=="tritonInelastic") || (process=="He3Inelastic")) return true; // wrong!!

	//ok for g4.10
	if((process =="ionInelastic") || (process =="alphaInelastic") || (process=="protonInelastic") || (process=="neutronInelastic") || (process=="dInelastic") || (process=="tInelastic") || (process=="He3Inelastic")) return true; 


	else return false;
	
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
G4VTrajectory* EventAction::GetParentTrajectory(G4int parentID)
{
  G4TrajectoryContainer* container = 
   
  G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
  if(container==0) return 0;
  size_t nTraj = container->size();
  for(size_t i=0;i<nTraj;i++)
  {
    G4VTrajectory* tr1 = (G4VTrajectory*)((*container)[i]);
    if(tr1->GetTrackID()==parentID) return tr1;
  }
  return 0;
}

*/
