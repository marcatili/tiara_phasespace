#include "HodoscopeSD.hh"
#include "HodoscopeHit.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"

HodoscopeSD::HodoscopeSD(G4String name)
:G4VSensitiveDetector(name)//, dhitsCollection(0), HCID(-1)
{
  G4String HCname;
  collectionName.insert(HCname="HodoscopehitsCollection");

  // For hits collection ID 
  HCID = -1;
}

HodoscopeSD::~HodoscopeSD(){;}

void HodoscopeSD::Initialize(G4HCofThisEvent* HCE)
{
   dhitsCollection = new HodoscopeHitsCollection(SensitiveDetectorName,collectionName[0]);

   if(HCID<0) 
	// { HCID = GetCollectionID(0); }
           { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(dhitsCollection); }
    HCE->AddHitsCollection(HCID,dhitsCollection);
}




G4bool HodoscopeSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  //Get energy deposited
  //G4double edep = aStep->GetTotalEnergyDeposit();
  //  if (edep==0.) return true;

  //get volume name
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());
  G4String PhysVolName = touchable->GetVolume()->GetName(); 
  //G4cout << "VOLUME NAME "<< PhysVolName <<G4endl;
  

  //G4int copyNo = touchable->GetVolume()->GetCopyNo();
   
 // G4ThreeVector worldPos = preStep->GetPosition();
	
  HodoscopeHit* newHit = new HodoscopeHit();
	//newHit->SetPixelNo( copyNo );	

        newHit->SetEnergy( preStep->GetKineticEnergy() );
        newHit->SetParticle( aStep->GetTrack()->GetDefinition() );
	newHit->SetTime (preStep->GetGlobalTime()); //time since the current event began
        newHit->SetPartPosition(preStep->GetPosition());
        newHit->SetPartDirection(preStep->GetMomentumDirection());

	newHit->SetVolumeName( PhysVolName );	
	newHit->SetTID( aStep->GetTrack()->GetTrackID() );
	newHit->SetPID( aStep->GetTrack()->GetParentID() );

        newHit->SetKinEnergy( aStep->GetTrack()->GetVertexKineticEnergy() );  
	
        if(aStep->GetTrack()->GetCreatorProcess()!=0)
    	 	newHit->SetProcess(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
        newHit->SetVertexVolumeName(aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName());

        //newHit->SetTrackVertexPos(aStep->GetTrack()->GetVertexPosition());
        //newHit->SetTrackVertexMomentum(aStep->GetTrack()->GetVertexMomentumDirection());

	//To select only particles that have just entered in the current volume, 
	//ie. is at the first step in the volume;   
	//the preStepPoint is at boundary (to select exiting particle, postStepPoint should be at boundary)

	//check what happens with and without if
	if(preStep->GetStepStatus() == fGeomBoundary){
		  dhitsCollection->insert( newHit );
	}
	//to stop tracking every particle (except primary protons) that has entered the phase space
	//this avoid scoring many times the same particles - and should make the simulation faster
	G4Track *aTrack=aStep->GetTrack();
	G4int particleid=aTrack->GetDefinition()->GetPDGEncoding();
  	if (particleid!=2212) aTrack->SetTrackStatus(fStopAndKill); 

  return true;
}

void HodoscopeSD::EndOfEvent(G4HCofThisEvent*)
{;}

void HodoscopeSD::clear()
{;} 

void HodoscopeSD::DrawAll()
{;} 

void HodoscopeSD::PrintAll()
{;} 
