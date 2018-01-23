#ifndef HodoscopeSD_h
#define HodoscopeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HodoscopeHit.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class HodoscopeSD : public G4VSensitiveDetector
{

  public:
      HodoscopeSD(G4String name);
      virtual ~HodoscopeSD();


      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      G4int HCID;

      HodoscopeHitsCollection *dhitsCollection;

  public:
};




#endif

