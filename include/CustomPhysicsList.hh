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
/// \file eventgenerator/particleGun/include/PhysicsList.hh
/// \brief Definition of the PhysicsList class
//
//
// $Id: PhysicsList.hh 68024 2013-03-13 13:42:01Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CustomPhysicsList_h
#define CustomPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
//class PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CustomPhysicsList: public G4VModularPhysicsList
{
  public:
    CustomPhysicsList();
   ~CustomPhysicsList();

  //protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess(); 
    virtual void SetCuts();   

  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void SetCutForProton(G4double);
        


  private:

  G4double fCutForGamma;
  G4double fCutForElectron;
  G4double fCutForPositron;
  G4double fCutForProton;

  G4VPhysicsConstructor*  fEmPhysicsList;
  G4VPhysicsConstructor*  fParticleList;
  G4VPhysicsConstructor*  fHadronPhys;
    
  G4VPhysicsConstructor*  fExtraPhys;
  G4VPhysicsConstructor*  fHadESPhys;
  G4VPhysicsConstructor*  fIonPhys;
  G4VPhysicsConstructor*  fStopPhys;
  G4VPhysicsConstructor*  fnTCutPhys;
  G4VPhysicsConstructor* fRadPhys;

  //PhysicsListMessenger* fMessenger;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif



