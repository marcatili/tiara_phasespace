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

//this physics list is based on  QGSP_BIC_EMY
//EMY stands for Emstandardphysics_option3
//
//QGSP 
//QGSP is the basic physics list applying the quark gluon string model for high energy interactions of protons, neutrons, pions, and Kaons and nuclei. The high energy interaction creates an exited nucleus, which is passed to the precompound model modeling the nuclear de-excitation.
//
//QGSP_BIC and QGSP_BIC_HP
//Like QGSP, but using Geant4 Binary cascade for primary protons and neutrons with energies below ~10GeV, thus replacing the use of the LEP model for protons and neutrons In comparison to the LEP model, Binary cascade better describes production of secondary particles produced in interactions of protons and neutrons with nuclei. Both lists, QGSP_BIC and QGSP_BIC_HP, also use the binary light ion cascade to model inelastic interaction of ions up to few GeV/nucleon with matter. The list QGSP_BIC_HP is like QGSP_BIC with the addition to use the data driven high precision neutron package (NeutronHP) to transport neutrons below 20 MeV down to thermal energies. 

//next: look at HP!!! to see if it makes some difference (neutron should be eliminated by tof, anyway)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CustomPhysicsList.hh"
//#include "PhysicsListMessenger.hh"


#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmProcessOptions.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ProductionCuts.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4IonPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CustomPhysicsList::CustomPhysicsList()
: G4VModularPhysicsList(),
fEmPhysicsList(0), fParticleList(0)//, fMessenger(0)
{ 

  G4bool HP_physics=true;	

  //defaultCutValue = 0.001*mm;
  defaultCutValue = 1.*mm;
  fCutForGamma     = defaultCutValue;
  fCutForElectron  = defaultCutValue;
  fCutForPositron  = defaultCutValue;
  fCutForProton    = defaultCutValue;
  verboseLevel    = 1;

  //fMessenger = new PhysicsListMessenger(this);

  // Particles
  fParticleList = new G4DecayPhysics();

  //EM physics
  fRadPhys = new G4RadioactiveDecayPhysics();

  //fEmPhysicsList = new G4EmStandardPhysics();
  fEmPhysicsList = new G4EmStandardPhysics_option3();

  // Synchroton Radiation & GN Physics
  fExtraPhys = new G4EmExtraPhysics();

  //hadron physics
   if (HP_physics) fHadronPhys = new G4HadronPhysicsQGSP_BIC_HP();
   else fHadronPhys = new G4HadronPhysicsQGSP_BIC();

  //Hadron Elastic scattering
   if (HP_physics)  fHadESPhys = new G4HadronElasticPhysicsHP( );
   else fHadESPhys = new G4HadronElasticPhysics( );

   // Ion Physics 
   if (HP_physics) fIonPhys = new G4IonPhysics();
   else fIonPhys = new G4IonBinaryCascadePhysics();   

    // Stopping Physics
    fStopPhys = new G4StoppingPhysics( );
    
    //Neutron tracking cut
    fnTCutPhys = new G4NeutronTrackingCut();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CustomPhysicsList::~CustomPhysicsList()
{ 
  //delete fMessenger;
  delete fParticleList;
  delete fEmPhysicsList;
  delete fHadronPhys;
  delete fExtraPhys;
  delete fHadESPhys;
  delete fIonPhys;
  delete fStopPhys;
  delete fnTCutPhys;
  delete fRadPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CustomPhysicsList::ConstructParticle()
{

    fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CustomPhysicsList::ConstructProcess()
{
  AddTransportation();

  fEmPhysicsList->ConstructProcess();
  fParticleList->ConstructProcess();
  fRadPhys->ConstructProcess();
  fHadronPhys->ConstructProcess();
  fExtraPhys->ConstructProcess();
  fHadESPhys->ConstructProcess();
  fIonPhys->ConstructProcess();
  fStopPhys->ConstructProcess();
  fnTCutPhys->ConstructProcess();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CustomPhysicsList::SetCuts()
{

  G4double HighPrecisionCutValue=1.*um;  
  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250.*eV, 1.*GeV);

  SetCutsWithDefault(); 
  //SetCutValue(fCutForGamma, "gamma");
  //SetCutValue(fCutForElectron, "e-");
  //SetCutValue(fCutForPositron, "e+");
  //SetCutValue(fCutForProton, "proton");


  // Production thresholds for detector regions
  G4ProductionCuts* cuts = new G4ProductionCuts;
  
  cuts->SetProductionCut(defaultCutValue);     //same cut for e+,e- photons...  
  /*
  cuts->SetProductionCut(HighPrecisionCutValue, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(HighPrecisionCutValue, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(HighPrecisionCutValue, G4ProductionCuts::GetIndex("e+"));
  */




  G4String regName ="DefaultRegionForTheWorld";
  G4Region* regDDet = G4RegionStore::GetInstance()->GetRegion(regName);  
  regDDet->SetProductionCuts(cuts);
/*
  G4String regName ="DiaReg";
  G4Region* regDDet = G4RegionStore::GetInstance()->GetRegion(regName);  
  regDDet->SetProductionCuts(cuts);
 
  regName ="AbsReg";
  G4Region* regADet = G4RegionStore::GetInstance()->GetRegion(regName);  
  regADet->SetProductionCuts(cuts);
*/

 regName ="PhaseSpaceReg";
  G4Region* regPhSP = G4RegionStore::GetInstance()->GetRegion(regName);  
  regPhSP->SetProductionCuts(cuts);

 regName ="headRegion";
  G4Region* reghead = G4RegionStore::GetInstance()->GetRegion(regName);  
  reghead->SetProductionCuts(cuts);




  if (verboseLevel>0) { 
	DumpCutValuesTable(); 
 }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CustomPhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CustomPhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CustomPhysicsList::SetCutForPositron(G4double cut)
{
  fCutForPositron = cut;
  SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CustomPhysicsList::SetCutForProton(G4double cut)
{
  fCutForProton = cut;
  SetParticleCuts(fCutForProton, G4Proton::Proton());
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

