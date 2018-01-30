
#ifndef HodoscopeHit_h
#define HodoscopeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

class HodoscopeHit : public G4VHit
{
  public:

      HodoscopeHit();
      ~HodoscopeHit();
      G4int operator==(const HodoscopeHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      
      G4ThreeVector position;
      G4ThreeVector momentum;
      G4double      energy;
      G4double      time;
      G4ParticleDefinition* particle;
      G4String volumename;
      G4int	trackid;
      G4int	parentid;

      G4double     ekin;
      G4String processname;
      G4String vtxvolumename;
      G4int	pixelnumber;
      G4ThreeVector trackvertex;
      G4ThreeVector trackvertexmom;
      G4ThreeVector partpos;
      G4ThreeVector partdir;
      G4double      vtime;

  public:
      
      inline void SetPosition(G4ThreeVector pos) { position=pos; }
      inline G4ThreeVector GetPosition() { return position; }
      inline void SetMomentum(G4ThreeVector mom) { momentum = mom; }
      inline G4ThreeVector GetMomentum(){ return momentum; }


      inline void SetEnergy(G4double ene){ energy = ene; }
      inline G4double GetEnergy(){ return energy; }
      inline void SetParticle(G4ParticleDefinition* pdef){ particle = pdef; }
      inline G4ParticleDefinition* GetParticle(){ return particle; }
      inline void SetTime(G4double val) { time = val; }
      inline G4double GetTime() const { return time; }  

      inline void SetVolumeName(G4String val) { volumename=val; }
      inline G4String GetVolumeName(){ return volumename; }
 
      inline void SetTID(G4int val) { trackid=val; }
      inline G4int GetTID(){ return trackid; }
      inline void SetPID(G4int val) { parentid=val; }
      inline G4int GetPID(){ return parentid; }

      inline void SetKinEnergy(G4double val){ ekin = val; }
      inline G4double GetKinEnergy(){ return ekin; }
     
      inline void SetProcess(G4String val) { processname=val; }
      inline G4String GetProcess(){ return processname; }

      inline void SetVertexVolumeName(G4String val) { vtxvolumename=val; }
      inline G4String GetVertexVolumeName(){ return vtxvolumename; }

      inline void SetPixelNo(G4int val) { pixelnumber=val; }
      inline G4int GetPixelNo(){ return pixelnumber; }
  
      inline void SetPartPosition(G4ThreeVector val) {partpos =val; }
      inline G4ThreeVector GetPartPosition(){ return partpos; }

      inline void SetPartDirection(G4ThreeVector val) {partdir =val; }
      inline G4ThreeVector GetPartDirection(){ return partdir; }

      inline void SetTrackVertexPos(G4ThreeVector val) {trackvertex =val; }
      inline G4ThreeVector GetTrackVertexPos(){ return trackvertex; }

      inline void SetTrackVertexMomentum(G4ThreeVector val) {trackvertexmom =val; }
      inline G4ThreeVector GetTrackVertexMomentum(){ return trackvertexmom; }

      inline void SetVertexTime(G4double val) { vtime = val; }
      inline G4double GetVertexTime() const { return vtime; }  


};

typedef G4THitsCollection<HodoscopeHit> HodoscopeHitsCollection;

extern G4Allocator<HodoscopeHit> HodoscopeHitAllocator;

inline void* HodoscopeHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HodoscopeHitAllocator.MallocSingle();
  return aHit;
}

inline void HodoscopeHit::operator delete(void *aHit)
{
  HodoscopeHitAllocator.FreeSingle((HodoscopeHit*) aHit);
}

#endif
