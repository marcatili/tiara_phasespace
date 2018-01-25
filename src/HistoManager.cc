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

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("tiara")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  
  analysis->SetFileName(fFileName);
  analysis->SetVerboseLevel(1);
  analysis->SetActivation(true);     //enable inactivation of histos, nTuples
 
  
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  //
  
  analysis->SetHistoDirectoryName("histo");  
  analysis->SetFirstHistoId(1);
    
  G4int id = analysis->CreateH1("h1.1","Gamma Prompt Energy", nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("h1.2","Phase Space impact time", nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
/* 
  id = analysis->CreateH1("h1.3","Triggering Particle Energy at Vertex", nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
 
 id = analysis->CreateH1("h1.4","Trigger time in diamond", nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
*/ 
// nTuples
  //
 
  analysis->SetNtupleDirectoryName("ntuple");
  analysis->SetFirstNtupleId(1);       
  analysis->CreateNtuple("ps", "Phase space Tuple");
  analysis->CreateNtupleIColumn(1, "evtnb");          //column 0 
  analysis->CreateNtupleDColumn(1, "Time");//1          
  analysis->CreateNtupleDColumn(1, "Ekin");          // 
  analysis->CreateNtupleDColumn(1, "x"); //
  analysis->CreateNtupleDColumn(1, "y"); //
  analysis->CreateNtupleDColumn(1, "z"); //
  analysis->CreateNtupleDColumn(1, "dx"); //
  analysis->CreateNtupleDColumn(1, "dy"); //
  analysis->CreateNtupleDColumn(1, "dz"); //
/*  
//vertex info 
  analysis->CreateNtupleDColumn(1, "vEkin");          // 
  analysis->CreateNtupleDColumn(1, "vx"); //
  analysis->CreateNtupleDColumn(1, "vy"); //
  analysis->CreateNtupleDColumn(1, "vz"); //
  analysis->CreateNtupleDColumn(1, "vdx"); //
  analysis->CreateNtupleDColumn(1, "vdy"); //
  analysis->CreateNtupleDColumn(1, "vdz"); //
*/

  analysis->CreateNtupleIColumn(1, "PID");         // //particle ID
  analysis->CreateNtupleIColumn(1, "TID");         // //track ID
  analysis->CreateNtupleIColumn(1, "VertexVolTag");  //
  analysis->CreateNtupleIColumn(1, "ProcessTag");    //
  analysis->CreateNtupleIColumn(1, "IsPG");          //
  analysis->CreateNtupleIColumn(1, "ParentID");          //

//  analysis->CreateNtupleDColumn(1, "tvx"); //track vertex x for triggering particle
//  analysis->CreateNtupleDColumn(1, "tvy"); //track vertex y for triggering particle
//  analysis->CreateNtupleDColumn(1, "tvz"); //track vertex z for triggering particle

  analysis->FinishNtuple(1);
 // analysis->SetNtupleActivation(true);         

//generated gamma prompt

  

//  analysis->SetNtupleDirectoryName("pg");  
  analysis->CreateNtuple("pg", "Prompt gamma Tuple"); //at vertex

  analysis->CreateNtupleIColumn(2, "evtnb"); //
  analysis->CreateNtupleIColumn(2, "TID"); //

  analysis->CreateNtupleDColumn(2, "PGx"); //PG vertex x 
  analysis->CreateNtupleDColumn(2, "PGy"); //PG vertex
  analysis->CreateNtupleDColumn(2, "PGz"); //PG vertex

  analysis->CreateNtupleDColumn(2, "PGdx");  
  analysis->CreateNtupleDColumn(2, "PGdy"); 
  analysis->CreateNtupleDColumn(2, "PGdz"); 

  analysis->CreateNtupleDColumn(2, "PGEkin"); //Ekin at vertex
  analysis->CreateNtupleDColumn(2, "PGtime"); //
  analysis->CreateNtupleIColumn(2, "ParentID");          //

  analysis->FinishNtuple(2);
 // analysis->SetNtupleActivation(true);             



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
