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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B4c::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace B4c
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Ar");
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  nistManager->FindOrBuildMaterial("G4_Cu") ;


  // define Elements
  G4Element* elH  = nistManager->FindOrBuildElement(1);
  G4Element* elC  = nistManager->FindOrBuildElement(6);
  G4Element* elSi = nistManager->FindOrBuildElement(14);
  G4Element* elO  = nistManager->FindOrBuildElement(8);


  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
          
  // Print materials        
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
          
  // Ar/CO2 mixture
  G4Material* Argon = nistManager->FindOrBuildMaterial("G4_Ar");
  G4Material* CarbonDioxide = nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4double mixtureDensity = (Argon->GetDensity() * 70/100.0 + CarbonDioxide->GetDensity() * 30/100.0) ;
  G4Material *ArCO2 = new G4Material("Ar/CO2",mixtureDensity,2) ;
  ArCO2->AddMaterial(Argon, 0.7) ;
  ArCO2->AddMaterial(CarbonDioxide, 0.3) ;   
  
  //FR4//
  //Epoxy (for FR4 )
  G4int numel(0), natoms(0);
  G4double fractionMass(0.);
  
  density = 1.2*g/cm3;
  G4Material* Epoxy = new G4Material("Epoxy" , density, numel=2);
  Epoxy->AddElement(elH, natoms=2);
  Epoxy->AddElement(elC, natoms=2);
  //SiO2 (Quarz)
  G4Material* SiO2 =  new G4Material("SiO2",density= 2.200*g/cm3, numel=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);
  //FR4 (Glass + Epoxy)
  density = 1.86*g/cm3;
  G4Material* FR4 = new G4Material("FR4"  , density, numel=2);
  FR4->AddMaterial(Epoxy, fractionMass=0.472);
  FR4->AddMaterial(SiO2, fractionMass=0.528);
              
 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  
  auto calorSizeXY  =  fNofPixelsXY * cellSizeXY; //10.*cm; //
  auto layerThickness = absoThickness + DriftThickness + gapThickness + PCBThickness;
  auto calorThickness = fNofLayers * layerThickness;
  auto worldSizeXY = 1.2 * calorSizeXY;
  auto worldSizeZ  = 1.2 * calorThickness;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  //auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  // auto gapMaterial = G4Material::GetMaterial("liquidArgon");
  //auto absorberMaterial = G4Material::GetMaterial("G4_Fe");
  //auto gapMaterial = G4Material::GetMaterial("G4_Ar");
  auto gapMaterial = G4Material::GetMaterial("Ar/CO2");
  auto absorberMaterial = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto PCBMaterial = G4Material::GetMaterial("FR4");
  auto driftMaterial = G4Material::GetMaterial("G4_Cu");
  
  //auto absorberMaterial = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  //auto PCBMaterial = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  //auto driftMaterial = G4Material::GetMaterial("G4_STAINLESS-STEEL");
    

  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World",          // its name
                  0, 0, 0);

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Calorimeter
  //
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size

  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter",     // its name
                  0, 0, 0);

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

    fScoringVolume = calorLV;

    //
    // Row plane - Y Slice
    //
    auto rowS
      = new G4Box("Row",           // its name
                   calorSizeXY/2, cellSizeXY/2, calorThickness/2); //its size

    auto rowLV
      = new G4LogicalVolume(
                   rowS,           // its solid
                   defaultMaterial,  // its material
                   "Row");         // its name

    new G4PVReplica(
                   "Row",          // its name
                   rowLV,          // its logical volume
                   calorLV,          // its mother
                   kYAxis,           // axis of replication
                   fNofPixelsXY,        // number of replica
                   cellSizeXY);  // witdth of replica
    
    //
    // Columns  - X Slice
    //
    auto cellS
      = new G4Box("Cell",           // its name
                   cellSizeXY/2, cellSizeXY/2, calorThickness/2); //its size

    auto cellLV
      = new G4LogicalVolume(
                   cellS,           // its solid
                   defaultMaterial,  // its material
                   "Cell");         // its name

    new G4PVReplica(
                   "Cell",          // its name
                   cellLV,          // its logical volume
                   rowLV,          // its mother
                   kXAxis,           // axis of replication
                   fNofPixelsXY,        // number of replica
                   cellSizeXY);  // witdth of replica
 
  //
  // Layer
  //
  auto layerS
    = new G4Box("Layer",           // its name
                 cellSizeXY/2, cellSizeXY/2, layerThickness/2); //its size

  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 cellLV,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers,        // number of replica
                 layerThickness);  // witdth of replica

  //
  // Absorber
  //
  auto absorberS
    = new G4Box("Abso",            // its name
                cellSizeXY/2, cellSizeXY/2, absoThickness/2); // its size

  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness/2 + absoThickness/2), // its position
                 absorberLV,       // its logical volume
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
  //
  // Copper Drift
  //
    auto DriftBoard
    = new G4Box("Drift",             // its name
                cellSizeXY/2, cellSizeXY/2, DriftThickness/2); // its size

  auto DriftBoardLV
    = new G4LogicalVolume(
                 DriftBoard,             // its solid
                 driftMaterial,      // its material
                 "driftBoardLV");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness/2 + absoThickness + DriftThickness/2), // its position
                 DriftBoardLV,            // its logical volume
                 "Drift",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Gap
  //
  auto gapS
    = new G4Box("Gap",             // its name
                cellSizeXY/2, cellSizeXY/2, gapThickness/2); // its size

  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-layerThickness/2 + absoThickness + DriftThickness + gapThickness/2), // its position
                 gapLV,            // its logical volume
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
                 
   //
   // RO board with PCB
   //
   
   auto ROboard
     = new G4Box("ReadOut",	//its name
     		 cellSizeXY/2, cellSizeXY/2, PCBThickness/2); //its size
     		 
   auto ROboardLV
    = new G4LogicalVolume(
                 ROboard,             // its solid
                 PCBMaterial,         // its material
                 "ROBoardLV");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-layerThickness/2 + absoThickness + DriftThickness + gapThickness + PCBThickness/2), // its position
                 ROboardLV,            // its logical volume
                 "ReadOut",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
   
  //
  // print parameters
  //
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName()
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(false);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  auto absVisAtt= new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  absVisAtt->SetVisibility(true);
  absorberLV->SetVisAttributes(absVisAtt);
    
  auto gapVisAtt= new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  gapVisAtt->SetVisibility(true);
  gapLV->SetVisAttributes(gapVisAtt);

  auto PCBVisAtt= new G4VisAttributes(G4Color(1.0,1.0,0.0, 0.5));
  PCBVisAtt->SetVisibility(true);
  ROboardLV->SetVisAttributes(PCBVisAtt);
  
  auto DriftVisAtt= new G4VisAttributes(G4Colour(0.0, 2.0, 0.0));
  DriftVisAtt->SetVisibility(true);
  DriftBoardLV->SetVisAttributes(DriftVisAtt);

  //
  // Always return the physical World
  //
  
   // Example of User Limits
    //
    // Below is an example of how to set tracking constraints in a given
    // logical volume
    //
    // Sets a max step length in the tracker region, with G4StepLimiter
    G4double maxStep = 0.5*gapThickness;
    fStepLimit = new G4UserLimits(maxStep);
    //logicSheet->SetUserLimits(fStepLimit);
    
  
  
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  //G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  //
  // Sensitive detectors
  //
  auto absoSD
    = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);
  SetSensitiveDetector("driftBoardLV",absoSD);
  SetSensitiveDetector("ROBoardLV",absoSD);

  auto gapSD
    = new CalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);
  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

}
