#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
//#include "HistoManager.hh"
#include "G4VTrajectory.hh"

class RunAction;


class EventAction : public G4UserEventAction
{
  public:
  
 EventAction(RunAction* runaction);
// EventAction();
    ~EventAction();

	

  public:
    G4bool IsInelastic(G4String process);
    G4VTrajectory* GetParentTrajectory(G4int parentID);
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
    //G4int hitsCollID;
      //G4int diaHCID; //
      //G4int absHCID; //
     G4int phspHCID; //
  
     //G4double edep;	
    // G4double edepTotMonoD=0.;	
    // G4double timeMonoD;	
    // G4double TMPtimeMonoD;	

  
    // G4String TrigPartName;
    // G4int TrigPartPIDMonoD;

    // G4double   TrEdep=0.;
	
    //  HistoManager* fHistoManager ;
     RunAction* fRunAction;

};

#endif

    
