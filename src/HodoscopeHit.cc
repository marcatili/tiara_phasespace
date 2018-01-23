#include "HodoscopeHit.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


G4Allocator<HodoscopeHit> HodoscopeHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HodoscopeHit::HodoscopeHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HodoscopeHit::~HodoscopeHit()
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4int HodoscopeHit::operator==(const HodoscopeHit &right) const
{
  return (this==&right) ? 1 : 0;
}


void HodoscopeHit::Draw()
{;}

void HodoscopeHit::Print()
{
//G4cout << "  HIT CLASS " << energy/MeV << " (MeV)" << G4endl;

}
