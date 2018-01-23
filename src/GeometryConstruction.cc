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

#include "GeometryConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
//#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "HodoscopeSD.hh"
#include "G4Region.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GeometryConstruction::GeometryConstruction()
: G4VUserDetectorConstruction(),
  fworld_phys(0),
  //fDia_phys(0),
  //fAbs_phys(0),
  fStepLimit(NULL)
{

    theMessenger = new GeometryConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GeometryConstruction::~GeometryConstruction()
{
  delete fStepLimit;
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* GeometryConstruction::Construct()
{
  ////////////////////////////////////////////////////////////////////////////////
  // Define materials.
  ////////////////////////////////////////////////////////////////////////////////
  G4double a, z, density; 
  G4int ncomponents;
  G4double fractionmass;

  
  G4String name, symbol;
  //G4NistManager* man = G4NistManager::Instance();

 //ELEMENTS
  G4Element* fH = new G4Element("H", "H", z=1., a=1.01*g/mole);
  G4Element* fC = new G4Element("C", "C", z=6., a=12.01*g/mole);
  G4Element* fN = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  G4Element* fO = new G4Element("O"  , "O", z=8., a= 16.00*g/mole);
  G4Element* fAr = new G4Element("Argon",    "Ar", z=18, a=  39.9480*g/mole);

  G4Element* fNa = new G4Element("Sodium","Na",z = 11., a=22.99*g/mole);
  G4Element* fMg = new G4Element("Magnesium","Mg",z = 12.,a=24.305*g/mole);
  G4Element* fP = new G4Element("Phosphorus","P",z = 15.,a=30.974*g/mole);
  G4Element* fS = new G4Element("Sulfur","S",z = 16.,a=32.064*g/mole);
  G4Element* fCl = new G4Element("Chlorine","Cl",z = 17.,a=35.453*g/mole);
  G4Element* fK = new G4Element("Potassium","K",z = 19.,a=39.098*g/mole);
  G4Element* fCa = new G4Element("Calcium","Ca",z = 20.,a=40.08*g/mole);
  G4Element* fFe  = new G4Element("Iron","Fe",z = 26.,a=55.85*g/mole);
  G4Element* fZn = new G4Element("Zinc","Zn",z = 30.,a=65.38*g/mole);
  G4Element* fRb = new G4Element("Rb","Rb",z = 37.,a=85.47 *g/mole);
  G4Element* fSr = new G4Element("Sr","Sr",z = 38.,a=87.62 *g/mole);
  G4Element* fZr = new G4Element("Zr","Zr",z = 40.,a=91.22 *g/mole);
  G4Element* fI = new G4Element("Iodine","I",z = 53.,a=126.90*g/mole);


//  G4Element* fCu = new G4Element("Copper",   "Cu", z=29, a=  63.5460*g/mole);
//  G4Element* fW = new G4Element(name="Tungsten" ,symbol="W",  z=74., 183.84*g/mole); //19.3 g/cm3
  G4Element* fPb = new G4Element(name="Lead"    ,symbol="Pb", z=82., 207.20*g/mole);  //11.35 g/cm3

//MATERIALS

	G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  //Air	
	G4Material* Air = 
	new G4Material("Air", density= 1.205*mg/cm3, ncomponents=4,
                      kStateGas, 293.*kelvin, 1.*atmosphere);
	Air->AddElement(fC, fractionmass=0.000124);                      
	Air->AddElement(fN, fractionmass=0.755267);
	Air->AddElement(fO, fractionmass=0.231782);
	Air->AddElement(fAr,fractionmass=0.012827);
 
  //Water
	G4Material* H2O = new G4Material(name="H2O", density= 1.000*g/cm3, ncomponents=2);
	H2O->AddElement(fH, ncomponents=2);
	H2O->AddElement(fO, ncomponents=1);
/*
  //DIAMOND
  	G4Material* Diamond = new G4Material("Diamond", density= 3.52*g/cm3, ncomponents=1);
  	Diamond->AddElement(fC, fractionmass=1.0);

  //Tungsten
	G4Material* Tungsten = new G4Material("Tungsten", density= 19.3*g/cm3, ncomponents=1);
	Tungsten->AddElement(fW, fractionmass=1.0);
  //Lead
	G4Material* Lead = new G4Material("Lead", density= 11.35*g/cm3, ncomponents=1);
	Lead->AddElement(fPb, fractionmass=1.0);
*/	

/*
  // MIRD soft tissue

  G4Material* soft = new G4Material("soft_tissue",density=0.9869 *g/cm3,16);
  soft->AddElement(fH,0.1047);
  soft->AddElement(fC,0.2302);
  soft->AddElement(fN,0.0234);
  soft->AddElement(fO,0.6321);
  soft->AddElement(fNa,0.0013);
  soft->AddElement(fMg,0.00015);
  soft->AddElement(fP,0.0024);
  soft->AddElement(fS,0.0022);
  soft->AddElement(fCl,0.0014);
  soft->AddElement(fK,0.0021);
  soft->AddElement(fFe,0.000063);
  soft->AddElement(fZn,0.000032);
  soft->AddElement(fRb,0.0000057);
  soft->AddElement(fSr,0.00000034);
  soft->AddElement(fZr,0.000008);
  soft->AddElement(fPb,0.00000016);
 */

 // MIRD Brain
    G4Material* matBrain = new G4Material("Brain",density=1.050*g/cm3,13);
    matBrain->AddElement(fH,0.1070);
    matBrain->AddElement(fC,0.1440);
    matBrain->AddElement(fN,0.0220);
    matBrain->AddElement(fO,0.7130);
    matBrain->AddElement(fNa,0.0020);
    matBrain->AddElement(fMg,0.0000);
    matBrain->AddElement(fP,0.0040);
    matBrain->AddElement(fS,0.0020);
    matBrain->AddElement(fCl,0.0030);
    matBrain->AddElement(fK,0.0030);
    matBrain->AddElement(fCa,0.0000);
    matBrain->AddElement(fFe,0.0000);
    matBrain->AddElement(fI,0.0000);


  // ICRP Cranium Spongiosa
    
    G4Material* matCranSpong = new G4Material("CraniumSpongiosa",density=1.245*g/cm3,13);
    matCranSpong->AddElement(fH,0.0810);
    matCranSpong->AddElement(fC,0.3170);
    matCranSpong->AddElement(fN,0.0280);
    matCranSpong->AddElement(fO,0.4510);
    matCranSpong->AddElement(fNa,0.0020);
    matCranSpong->AddElement(fMg,0.0010);
    matCranSpong->AddElement(fP,0.0370);
    matCranSpong->AddElement(fS,0.0030);
    matCranSpong->AddElement(fCl,0.0010);
    matCranSpong->AddElement(fK,0.0010);
    matCranSpong->AddElement(fCa,0.0780);
    matCranSpong->AddElement(fFe,0.0000);
    matCranSpong->AddElement(fI,0.0000);

  ////////////////////////////////////////////////////////////////////////////////
  // END materials.
  ////////////////////////////////////////////////////////////////////////////////


   //------------------------------ 
   // Regions
   //-----------------------------

 //   G4Region* DiamondReg = new G4Region("DiaReg");
 //   G4Region* AbsorberReg = new G4Region("AbsReg");
      G4Region* PhSpace = new G4Region("PhaseSpaceReg");
      G4Region* headRegion = new G4Region("headRegion");

    //check overlaps
    G4bool pSurfChk=false;

 
  // Define size of world and volumes in it.

  G4double world_x = 50*cm;
  G4double world_y = 50*cm;
  G4double world_z = 50*cm;

	
 //skull 
  G4double skull_thickness = 7.*mm;
  G4double skull_Rmax = 10.*cm;
  G4double skull_Rmin = skull_Rmax-skull_thickness ;
 
	

 //head
  G4double head_R= skull_Rmin;

/*
  //Diamond Detector
  G4double Dia_x = (500./2.)*um; 
  G4double Dia_y = 5./2.*mm;
  G4double Dia_z = 5./2.*mm;

  
 //Absorber

  Abs_x = (800./2.)*um; //to be set by macro
  G4double Abs_y = Dia_y;
  G4double Abs_z = Dia_z;


  G4double Abs_posR= skull_Rmax+Abs_x; 
  G4double Dia_posR= skull_Rmax + 2.*Abs_x + Dia_x; 

  


  G4cout<<"_______GEOMETRY_______"<<G4endl;	
  G4cout<<"Abs_x = "<<Abs_x/um<<" um"<<G4endl;	
*/

  //WORLD
	  G4Box * world=  new G4Box("world", world_x, world_y, world_z);
	  G4LogicalVolume * world_log = new G4LogicalVolume(world, Air,"world_L",0,0,0);
	  fworld_phys=  new G4PVPlacement(0,G4ThreeVector(),"world_phys", world_log,0,false,0);


//////////////////////////////////////////////////////////////////////////////


   //head

         G4Sphere *head = new G4Sphere("head", 0., head_R,  0, 2*pi, 0, pi);
         head_log = new G4LogicalVolume(head,matBrain,"head_log",0,0,0);
	 head_phys = new G4PVPlacement(0,G4ThreeVector(0. ,0.,0.),"head_phys", head_log,fworld_phys,false,0, pSurfChk);

         G4Sphere *skull = new G4Sphere("skull", skull_Rmin, skull_Rmax,  0, 2*pi, 0, pi);
         skull_log = new G4LogicalVolume(skull,matCranSpong,"head_log",0,0,0);
	 skull_phys = new G4PVPlacement(0,G4ThreeVector(0. ,0.,0.),"skull_phys", skull_log,fworld_phys,false,0, pSurfChk);

	 // Region
	  head_log->SetRegion(headRegion);
	  headRegion->AddRootLogicalVolume(head_log);

 	  skull_log->SetRegion(headRegion);
	  headRegion->AddRootLogicalVolume(skull_log);

    //Detectors
/*

  	////////////////////////////////////
	//Diamond 		
	////////////////////////////////////
	 G4Box * Dia = new G4Box("Dia", Dia_x, Dia_y, Dia_z);	
         Dia_log = new G4LogicalVolume(Dia,Diamond,"Dia_log",0,0,0);

	 Dia_log->SetRegion(DiamondReg);
         DiamondReg -> AddRootLogicalVolume(Dia_log);

	 G4double theta=0.;
	 G4int cpn=0;
	 G4RotationMatrix *rotm;;
	 G4double AngSpacing=pi/6.;


	//phi=0 only one detector
	G4double phi=0.;
	G4double w = std::cos(phi);
	rotm  = new G4RotationMatrix();
	rotm->rotateY((pi/2.)-phi); //this is the good one
		
       	 new G4PVPlacement(rotm, G4ThreeVector(0, 0,Dia_posR*std::cos(phi)), "Dia_phys", Dia_log,fworld_phys,false,cpn, pSurfChk);  
				
	G4cout<<"Copy Nb "<<cpn<<"  Detector position: "<< G4ThreeVector( 0,0,Dia_posR*std::cos(phi)) << G4endl;
	cpn++;	



   	 for(phi=pi/6.; phi<=pi/2.; phi=phi+AngSpacing){
	 	for(theta=0.; theta<2.*pi; theta=theta+AngSpacing){ 

			//if (std::abs((phi-(pi/2.))<0.1) && (std::abs(theta - pi)<0.1)){} else { 
			if ((phi==pi/2.) && (std::abs(theta - pi)<0.1)){} else { 
	 			w = std::cos(phi);
				rotm  = new G4RotationMatrix(G4ThreeVector(0, 0, w), -theta);//this is the good one
				rotm->rotateY((pi/2.)-phi); //this is the good one
		
		        	 new G4PVPlacement(rotm, G4ThreeVector( Dia_posR*std::cos(theta)*std::sin(phi),	
					Dia_posR*std::sin(theta)*std::sin(phi),Dia_posR*std::cos(phi)), 
					"Dia_phys", Dia_log,fworld_phys,false,cpn, pSurfChk);  
				
				G4cout<<"Copy Nb "<<cpn<<"  Detector position: "<< G4ThreeVector( Dia_posR*std::cos(theta)*std::sin(phi),	
					Dia_posR*std::sin(theta)*std::sin(phi),Dia_posR*std::cos(phi)) << G4endl;

				cpn++;	
			}
		}
	}

	for(phi=pi/2.+(pi/6.); phi<(pi-0.1); phi=phi+AngSpacing){
		for(theta=0.; theta<2.*pi; theta=theta+AngSpacing){ 

 			w = std::cos(phi);
			rotm  = new G4RotationMatrix(G4ThreeVector(0, 0, w), theta);//this is the good one
			rotm->rotateY((pi/2.)-phi); //this is the good one

		         new G4PVPlacement(rotm, G4ThreeVector( Dia_posR*std::cos(theta)*std::sin(phi),
				Dia_posR*std::sin(theta)*std::sin(phi),Dia_posR*std::cos(phi)), 
				"Dia_phys", Dia_log,fworld_phys,false,cpn, pSurfChk);  

				G4cout<<"Copy Nb "<<cpn<<"  Detector position: "<< G4ThreeVector( Dia_posR*std::cos(theta)*std::sin(phi),	
					Dia_posR*std::sin(theta)*std::sin(phi),Dia_posR*std::cos(phi))<<G4endl;


	  		 cpn++;	
		}
	}

	//phi=pi only one detector
	phi=pi;
	w = std::cos(phi);
	rotm  = new G4RotationMatrix();
	rotm->rotateY((pi/2.)-phi); //this is the good one
		
       	 new G4PVPlacement(rotm, G4ThreeVector(0, 0, Dia_posR*std::cos(phi)), "Dia_phys", Dia_log,fworld_phys,false,cpn, pSurfChk);  
				
	G4cout<<"Copy Nb "<<cpn<<"  Detector position: "<< G4ThreeVector( 0,0,Dia_posR*std::cos(phi)) << G4endl;




	G4cout<<"Number of pixel defined = "<<cpn<<G4endl;
   	G4cout<<"_____________________"<<G4endl;	


	////////////////////////////////////
        //Absorbers
	////////////////////////////////////

	 G4Box *Abs = new G4Box("Abs", Abs_x, Abs_y, Abs_z);	
         Abs_log = new G4LogicalVolume(Abs, Tungsten,"Abs_log", 0, 0, 0);
  	
	 Abs_log->SetRegion(AbsorberReg);
         AbsorberReg -> AddRootLogicalVolume(Abs_log);

        // fAbs_phys = new G4PVPlacement(0,G4ThreeVector(Abs_posR ,0.,0.),"Abs_phys", Abs_log,fworld_phys,false,0, pSurfChk);
	

	cpn=0;
	phi=0.;
	theta=0.;


	//phi=0 only one detector
	phi=0.;
	w = std::cos(phi);
	rotm  = new G4RotationMatrix();
	rotm->rotateY((pi/2.)-phi); //this is the good one
		
       	 new G4PVPlacement(rotm, G4ThreeVector(0, 0,Abs_posR*std::cos(phi)), "Abs_phys", Abs_log,fworld_phys,false,cpn, pSurfChk); 



  	for(phi=pi/6.; phi<=pi/2.; phi=phi+AngSpacing){
	 	for(theta=0.; theta<2.*pi; theta=theta+AngSpacing){ 

			//if (std::abs((phi-(pi/2.))<0.1) && (std::abs(theta - pi)<0.1)){} else { 
			if ((phi==pi/2.) && (std::abs(theta - pi)<0.1)){} else { 
	 			w = std::cos(phi);
				rotm  = new G4RotationMatrix(G4ThreeVector(0, 0, w), -theta);//this is the good one
				rotm->rotateY((pi/2.)-phi); //this is the good one
		
		        	 new G4PVPlacement(rotm, G4ThreeVector( Abs_posR*std::cos(theta)*std::sin(phi),	
					Abs_posR*std::sin(theta)*std::sin(phi),Abs_posR*std::cos(phi)), 
					"Abs_phys", Abs_log,fworld_phys,false,cpn, pSurfChk);  
				cpn++;	
			}
		}
	}

	for(phi=pi/2.+(pi/6.); phi<(pi-0.1); phi=phi+AngSpacing){
		for(theta=0.; theta<2.*pi; theta=theta+AngSpacing){ 

 			w = std::cos(phi);
			rotm  = new G4RotationMatrix(G4ThreeVector(0, 0, w), theta);//this is the good one
			rotm->rotateY((pi/2.)-phi); //this is the good one

		         new G4PVPlacement(rotm, G4ThreeVector( Abs_posR*std::cos(theta)*std::sin(phi),
				Abs_posR*std::sin(theta)*std::sin(phi),Abs_posR*std::cos(phi)), 
				"Abs_phys", Abs_log,fworld_phys,false,cpn, pSurfChk); 
	  		 cpn++;	
		}
	}

	//phi=pi only one detector
	phi=pi;
	w = std::cos(phi);
	rotm  = new G4RotationMatrix();
	rotm->rotateY((pi/2.)-phi); //this is the good one
		
       	 new G4PVPlacement(rotm, G4ThreeVector(0, 0,Abs_posR*std::cos(phi)), "Abs_phys", Abs_log,fworld_phys,false,cpn, pSurfChk); 
*/



//////////////////
// PHASE SPACE
//////////////////

  //G4double 
  //ShellThickness = 5.*mm;
  G4cout<<"Phase space thickness   "<< ShellThickness/mm<< " mm "<<G4endl;		
  G4double PhSp_Rmin=skull_Rmax;
  G4double PhSp_Rmax=PhSp_Rmin+ShellThickness;


  G4Sphere *PhSpShell = new G4Sphere("PhSpShell", PhSp_Rmin, PhSp_Rmax,  0, 2*pi, 0, pi);
         PhSpShell_log = new G4LogicalVolume(PhSpShell,Vacuum,"PhSpShell_log",0,0,0);
	 PhSpShell_phys = new G4PVPlacement(0,G4ThreeVector(0. ,0.,0.),"PhSpShell_phys", PhSpShell_log,fworld_phys,false,0, pSurfChk);

  
  // Region
  PhSpShell_log->SetRegion(PhSpace);
  PhSpace->AddRootLogicalVolume(PhSpShell_log);






//--------- Visualization attributes -------------------------------
  world_log->SetVisAttributes(G4VisAttributes::Invisible);
 
  G4VisAttributes* aVisAtt= new G4VisAttributes(G4Colour(0,0,1));
  //Abs_log->SetVisAttributes(aVisAtt);
  //Abs_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  //G4VisAttributes* bVisAtt= new G4VisAttributes(G4Colour(1.0,.0,1.0));
  //bVisAtt->SetForceSolid(true);
  //Dia_log->SetVisAttributes(bVisAtt);
  //Dia_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  //G4VisAttributes* cVisAtt= new G4VisAttributes(G4Colour(1.0,1.,1.0));
  //cVisAtt->SetForceSolid(true);
  head_log->SetVisAttributes(G4VisAttributes::Invisible);
  //head_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  skull_log->SetVisAttributes(aVisAtt);

  G4VisAttributes* phsp_vis = new G4VisAttributes(1);
  G4Colour mycolor(G4Colour::Yellow());
  phsp_vis->SetColor(mycolor);
  PhSpShell_log->SetVisAttributes(phsp_vis);	

// set tracking constraints in a given logical volume

// Sets a max step length 

   G4double  maxStep = 2.*ShellThickness;
   fStepLimit = new G4UserLimits(maxStep);
   PhSpShell_log->SetUserLimits(fStepLimit);


  /*
   G4double  maxStep = 0.1*Dia_x;
   fStepLimit = new G4UserLimits(maxStep);
   Dia_log->SetUserLimits(fStepLimit);
   Abs_log->SetUserLimits(fStepLimit);

   /// Set additional contraints on the track, with G4UserSpecialCuts
   ///
   /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
   /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
   ///                                           maxLength,
   ///                                           maxTime,
   ///                                           minEkin));
 
*/




  return fworld_phys;
}



void GeometryConstruction::ConstructSDandField()
{
    // sensitive detectors -----------------------------------------------------
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;

    G4VSensitiveDetector* PhaseSpaceSD = new HodoscopeSD(SDname="/phasespace");
    SDman->AddNewDetector(PhaseSpaceSD);
    PhSpShell_log->SetSensitiveDetector(PhaseSpaceSD);

/*
    G4VSensitiveDetector* DiamondSD = new HodoscopeSD(SDname="/diamond");
    SDman->AddNewDetector(DiamondSD);
    Dia_log->SetSensitiveDetector(DiamondSD);


    G4VSensitiveDetector* AbsorberSD = new HodoscopeSD(SDname="/absorber");
    SDman->AddNewDetector(AbsorberSD);
    Abs_log->SetSensitiveDetector(AbsorberSD);
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void GeometryConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
