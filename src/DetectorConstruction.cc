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
  nistManager->FindOrBuildMaterial("G4_AIR") ;


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
  G4double mixtureDensity = (Argon->GetDensity() * 93/100.0 + CarbonDioxide->GetDensity() * 7/100.0) ;
  G4Material *ArCO2 = new G4Material("Ar/CO2",mixtureDensity,2) ;
  ArCO2->AddMaterial(Argon, 0.93) ;
  ArCO2->AddMaterial(CarbonDioxide, 0.07) ;   
  
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
  
  auto calorSizeXY1  =  fNofPixelsXY1 * cellSizeXY; //10.*cm; //
  auto calorSizeXY2  =  fNofPixelsXY2 * cellSizeXY;
  
  auto layerThickness1 = airGapThickness + absoThickness1 + airGapThickness + DriftThickness + gapThickness + PCBThickness;
  auto calorThickness1 = fNofLayers1 * layerThickness1;
  
  auto layerThickness2 = airGapThickness + absoThickness2 + airGapThickness + DriftThickness + gapThickness + PCBThickness;
  auto calorThickness2 = fNofLayers2 * layerThickness2;
  
  auto worldSizeXY = 1.9 * std::max(calorSizeXY1, calorSizeXY2);
  auto worldSizeZ  = 1.9 * (calorThickness1 + calorThickness2);
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
  auto air = G4Material::GetMaterial("G4_AIR");
  
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
  auto calorimeterS1
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY1/2, calorSizeXY1/2, calorThickness1/2); // its size

  auto calorLV1
    = new G4LogicalVolume(
                 calorimeterS1,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter1",     // its name
                  0, 0, 0);

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV1,          // its logical volume
                 "Calorimeter1",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

    fScoringVolume = calorLV1;

    //
    // Row plane - Y Slice
    //
    auto rowS1
      = new G4Box("Row1",           // its name
                   calorSizeXY1/2, cellSizeXY/2, calorThickness1/2); //its size

    auto rowLV1
      = new G4LogicalVolume(
                   rowS1,           // its solid
                   defaultMaterial,  // its material
                   "Row1");         // its name

    new G4PVReplica(
                   "Row1",          // its name
                   rowLV1,          // its logical volume
                   calorLV1,          // its mother
                   kYAxis,           // axis of replication
                   fNofPixelsXY1,        // number of replica
                   cellSizeXY);  // witdth of replica
    
    //
    // Columns  - X Slice
    //
    auto cellS1
      = new G4Box("Cell1",           // its name
                   cellSizeXY/2, cellSizeXY/2, calorThickness1/2); //its size

    auto cellLV1
      = new G4LogicalVolume(
                   cellS1,           // its solid
                   defaultMaterial,  // its material
                   "Cell1");         // its name

    new G4PVReplica(
                   "Cell1",          // its name
                   cellLV1,          // its logical volume
                   rowLV1,          // its mother
                   kXAxis,           // axis of replication
                   fNofPixelsXY1,        // number of replica
                   cellSizeXY);  // witdth of replica
 
  //
  // Layer
  //
  auto layerS1
    = new G4Box("Layer",           // its name
                 cellSizeXY/2, cellSizeXY/2, layerThickness1/2); //its size

  auto layerLV1
    = new G4LogicalVolume(
                 layerS1,           // its solid
                 defaultMaterial,  // its material
                 "Layer1");         // its name

  new G4PVReplica(
                 "Layer1",          // its name
                 layerLV1,          // its logical volume
                 cellLV1,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers1,        // number of replica
                 layerThickness1);  // witdth of replica


   //
  // Air gap
  //
  auto airGap1
    = new G4Box("airGap1",            // its name
                cellSizeXY/2, cellSizeXY/2, airGapThickness/2); // its size

  auto airGapLV1
    = new G4LogicalVolume(
                 airGap1,        // its solid
                 air,		 // its material
                 "airGapLV1");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness1/2 + airGapThickness/2), // its position
                 airGapLV1,       // its logical volume
                 "airGapLV1",           // its name
                 layerLV1,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Absorber
  //
  auto absorberS1
    = new G4Box("Abso1",            // its name
                cellSizeXY/2, cellSizeXY/2, absoThickness1/2); // its size

  auto absorberLV1
    = new G4LogicalVolume(
                 absorberS1,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV1");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness1/2 + airGapThickness + absoThickness1/2), // its position
                 absorberLV1,       // its logical volume
                 "Abso1",           // its name
                 layerLV1,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
                 
      //
  // Air gap
  //
  auto airGap1_2
    = new G4Box("airGap1_2",            // its name
                cellSizeXY/2, cellSizeXY/2, airGapThickness/2); // its size

  auto airGapLV1_2
    = new G4LogicalVolume(
                 airGap1_2,        // its solid
                 air,		 // its material
                 "airGapLV1_2");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness1/2 + airGapThickness + absoThickness1 + airGapThickness/2), // its position
                 airGapLV1_2,       // its logical volume
                 "airGapLV1_2",           // its name
                 layerLV1,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
                 
  //
  // Drift (3 mm FR4 + 13 um Copper)
  //
    auto DriftBoard1
    = new G4Box("Drift1",             // its name
                cellSizeXY/2, cellSizeXY/2, DriftThickness/2); // its size

  auto DriftBoardLV1
    = new G4LogicalVolume(
                 DriftBoard1,             // its solid
                 PCBMaterial,      // its material
                 "driftBoardLV1");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness1/2 + airGapThickness + absoThickness1 + airGapThickness + DriftThickness/2), // its position
                 DriftBoardLV1,            // its logical volume
                 "Drift1",            // its name
                 layerLV1,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Gap
  //
  auto gapS1
    = new G4Box("Gap",             // its name
                cellSizeXY/2, cellSizeXY/2, gapThickness/2); // its size

  auto gapLV1
    = new G4LogicalVolume(
                 gapS1,             // its solid
                 gapMaterial,      // its material
                 "GapLV1");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-layerThickness1/2 + airGapThickness + absoThickness1 + airGapThickness + DriftThickness + gapThickness/2), // its position
                 gapLV1,            // its logical volume
                 "Gap1",            // its name
                 layerLV1,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
                 
   //
   // RO board with PCB
   //
   
   auto ROboard1
     = new G4Box("ReadOut1",	//its name
     		 cellSizeXY/2, cellSizeXY/2, PCBThickness/2); //its size
     		 
   auto ROboardLV1
    = new G4LogicalVolume(
                 ROboard1,             // its solid
                 PCBMaterial,         // its material
                 "ROBoardLV1");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-layerThickness1/2 + airGapThickness + absoThickness1 + airGapThickness + DriftThickness + gapThickness + PCBThickness/2), // its position
                 ROboardLV1,            // its logical volume
                 "ReadOut1",            // its name
                 layerLV1,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
                 
                 
  /////////////
  //different absorber thickness 
  /////////////
  
  //
  // Calorimeter
  //
  auto calorimeterS2
    = new G4Box("Calorimeter2",     // its name
                 calorSizeXY2/2, calorSizeXY2/2, calorThickness2/2); // its size

  auto calorLV2
    = new G4LogicalVolume(
                 calorimeterS2,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter2",     // its name
                  0, 0, 0);

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,0, calorThickness1/2 + calorThickness2/2),  // at (0,0,0)
                 calorLV2,          // its logical volume
                 "Calorimeter2",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

    //fScoringVolume = calorLV;

    //
    // Row plane - Y Slice
    //
    auto rowS2
      = new G4Box("Row2",           // its name
                   calorSizeXY2/2, cellSizeXY/2, calorThickness2/2); //its size

    auto rowLV2
      = new G4LogicalVolume(
                   rowS2,           // its solid
                   defaultMaterial,  // its material
                   "Row2");         // its name

    new G4PVReplica(
                   "Row2",          // its name
                   rowLV2,          // its logical volume
                   calorLV2,          // its mother
                   kYAxis,           // axis of replication
                   fNofPixelsXY2,        // number of replica
                   cellSizeXY);  // witdth of replica
    
    //
    // Columns  - X Slice
    //
    auto cellS2
      = new G4Box("Cell2",           // its name
                   cellSizeXY/2, cellSizeXY/2, calorThickness2/2); //its size

    auto cellLV2
      = new G4LogicalVolume(
                   cellS2,           // its solid
                   defaultMaterial,  // its material
                   "Cell2");         // its name

    new G4PVReplica(
                   "Cell2",          // its name
                   cellLV2,          // its logical volume
                   rowLV2,          // its mother
                   kXAxis,           // axis of replication
                   fNofPixelsXY2,        // number of replica
                   cellSizeXY);  // witdth of replica
 
  //
  // Layer
  //
  auto layerS2
    = new G4Box("Layer2",           // its name
                 cellSizeXY/2, cellSizeXY/2, layerThickness2/2); //its size

  auto layerLV2
    = new G4LogicalVolume(
                 layerS2,           // its solid
                 defaultMaterial,  // its material
                 "Layer2");         // its name

  new G4PVReplica(
                 "Layer2",          // its name
                 layerLV2,          // its logical volume
                 cellLV2,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers2,        // number of replica
                 layerThickness2);  // witdth of replica
                 
     //
  // Air gap
  //
  auto airGap2
    = new G4Box("airGap2",            // its name
                cellSizeXY/2, cellSizeXY/2, airGapThickness/2); // its size

  auto airGapLV2
    = new G4LogicalVolume(
                 airGap2,        // its solid
                 air,		 // its material
                 "airGapLV2");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness2/2 + airGapThickness/2), // its position
                 airGapLV2,       // its logical volume
                 "airGapLV2",           // its name
                 layerLV2,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


  //
  // Absorber
  //
  auto absorberS2
    = new G4Box("Abso2",            // its name
                cellSizeXY/2, cellSizeXY/2, absoThickness2/2); // its size

  auto absorberLV2
    = new G4LogicalVolume(
                 absorberS2,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV2");        // its name

//-layerThickness1/2 + absoThickness1 + DriftThickness + gapThickness + PCBThickness/2
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness2/2 + airGapThickness + absoThickness2/2), // its position
                 //G4ThreeVector(0., 0., -layerThickness1/2 + absoThickness1 + DriftThickness + gapThickness + PCBThickness -layerThickness2/2 + absoThickness2/2), // its position
                 absorberLV2,       // its logical volume
                 "Abso2",           // its name
                 layerLV2,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
  
      //
  // Air gap
  //
  auto airGap2_2
    = new G4Box("airGap2_2",            // its name
                cellSizeXY/2, cellSizeXY/2, airGapThickness/2); // its size

  auto airGapLV2_2
    = new G4LogicalVolume(
                 airGap2_2,        // its solid
                 air,		 // its material
                 "airGapLV2_2");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness2/2 + airGapThickness + absoThickness2 + airGapThickness/2), // its position
                 airGapLV2_2,       // its logical volume
                 "airGapLV2_2",           // its name
                 layerLV2,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
  
  //
  // Copper Drift
  //
    auto DriftBoard2
    = new G4Box("Drift2",             // its name
                cellSizeXY/2, cellSizeXY/2, DriftThickness/2); // its size

  auto DriftBoardLV2
    = new G4LogicalVolume(
                 DriftBoard2,             // its solid
                 PCBMaterial,      // its material
                 "driftBoardLV2");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness2/2 + airGapThickness + absoThickness2 + airGapThickness + DriftThickness/2), // its position
                 //G4ThreeVector(0., 0., -layerThickness1/2 + absoThickness1 + DriftThickness + gapThickness + PCBThickness -layerThickness2/2 + absoThickness2 + DriftThickness/2), // its position
                 DriftBoardLV2,            // its logical volume
                 "Drift2",            // its name
                 layerLV2,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Gap
  //
  auto gapS2
    = new G4Box("Gap2",             // its name
                cellSizeXY/2, cellSizeXY/2, gapThickness/2); // its size

  auto gapLV2
    = new G4LogicalVolume(
                 gapS2,             // its solid
                 gapMaterial,      // its material
                 "GapLV2");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness2/2 + airGapThickness + absoThickness2 + airGapThickness + DriftThickness + gapThickness/2), // its position
                 //G4ThreeVector(0., 0.,-layerThickness1/2 + absoThickness1 + DriftThickness + gapThickness + PCBThickness -layerThickness2/2 + absoThickness2  + DriftThickness + gapThickness/2), // its position
                 gapLV2,            // its logical volume
                 "Gap2",            // its name
                 layerLV2,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
                 
   //
   // RO board with PCB
   //
   
   auto ROboard2
     = new G4Box("ReadOut2",	//its name
     		 cellSizeXY/2, cellSizeXY/2, PCBThickness/2); //its size
     		 
   auto ROboardLV2
    = new G4LogicalVolume(
                 ROboard2,             // its solid
                 PCBMaterial,         // its material
                 "ROBoardLV2");        // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -layerThickness2/2 + airGapThickness + absoThickness2 + airGapThickness + DriftThickness + gapThickness + PCBThickness/2), // its position
                // G4ThreeVector(0., 0.,-layerThickness1/2 + absoThickness1 + DriftThickness + gapThickness + PCBThickness -layerThickness2/2 + absoThickness2 + DriftThickness + gapThickness + PCBThickness/2), // its position
                 ROboardLV2,            // its logical volume
                 "ReadOut2",            // its name
                 layerLV2,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
                 
   
  //
  // print parameters
  //
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers1 << " layers of: [ "
    << absoThickness1/mm << "mm of " << absorberMaterial->GetName()
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(false);
  calorLV1->SetVisAttributes(simpleBoxVisAtt);
  calorLV2->SetVisAttributes(simpleBoxVisAtt);

  auto absVisAtt= new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  absVisAtt->SetVisibility(true);
  absorberLV1->SetVisAttributes(absVisAtt);
  absorberLV2->SetVisAttributes(absVisAtt);
    
  auto gapVisAtt= new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  gapVisAtt->SetVisibility(true);
  gapLV1->SetVisAttributes(gapVisAtt);
  gapLV2->SetVisAttributes(gapVisAtt);

  auto PCBVisAtt= new G4VisAttributes(G4Color(1.0,1.0,0.0, 0.5));
  PCBVisAtt->SetVisibility(true);
  ROboardLV1->SetVisAttributes(PCBVisAtt);
  ROboardLV2->SetVisAttributes(PCBVisAtt);
  
  auto DriftVisAtt= new G4VisAttributes(G4Colour(0.0, 2.0, 0.0));
  DriftVisAtt->SetVisibility(true);
  DriftBoardLV1->SetVisAttributes(DriftVisAtt);
  DriftBoardLV2->SetVisAttributes(DriftVisAtt);

  //
  // Always return the physical World
  //
  
   // Example of User Limits
    //
    // Below is an example of how to set tracking constraints in a given
    // logical volume
    //
    // Sets a max step length in the tracker region, with G4StepLimiter
 //   G4double maxStep = 0.5*gapThickness;
 //   fStepLimit= new G4UserLimits(maxStep);
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
    = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection");
   // = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers1);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV1",absoSD);
  SetSensitiveDetector("driftBoardLV1",absoSD);
  SetSensitiveDetector("ROBoardLV1",absoSD);
  
  SetSensitiveDetector("AbsoLV2",absoSD);
  SetSensitiveDetector("driftBoardLV2",absoSD);
  SetSensitiveDetector("ROBoardLV2",absoSD);

  auto gapSD
      = new CalorimeterSD("GapSD", "GapHitsCollection");
      //= new CalorimeterSD("GapSD", "GapHitsCollection", fNofLayers1);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV1",gapSD);
  SetSensitiveDetector("GapLV2",gapSD);
  //SetSensitiveDetector("GapLV1",gapSD);
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
/*
void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}
*/
}
