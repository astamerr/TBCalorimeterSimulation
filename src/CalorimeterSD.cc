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
/// \file CalorimeterSD.cc
/// \brief Implementation of the B4c::CalorimeterSD class

#include "CalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace B4c
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(const G4String& name,
                             const G4String& hitsCollectionName,
                             G4int nofLayers)
 : G4VSensitiveDetector(name),
   fNofLayers(nofLayers)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection
    = new CalorHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

  // Create hits
  // fill calorimeter hits with zero energy deposition
  // nPixelX*nPixelY*(fNofLayers for layers + one more for total sums)
 for (auto column=0; column<fNofPixelsXY; column++) {
     for (auto row=0; row<fNofPixelsXY; row++) {
         for (G4int i=0; i<fNofLayers; i++ ) {
             fHitsCollection->insert(new CalorHit());
         }
     }
  }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CalorimeterSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{

    // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
 
  
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  //G4cout<< edep<< G4endl;
  analysisManager->FillH1(1, edep);


  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
   
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
  G4int particlePDG = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  
  
  G4int stepID = step->GetTrack()->GetCurrentStepNumber();
  G4int trackID = step->GetTrack()->GetTrackID() ;			// track ID
  G4double edepStep = step->GetTotalEnergyDeposit();
  G4ThreeVector Position = step->GetPreStepPoint()->GetPosition();
  
  if ( edep==0. && stepLength == 0. ) return false; 

  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  // Get calorimeter cell id
    //placement tree: what’s called the touchable’s “history”
    //touchable->GetCopyNumber(), you can get the cell copy number ?
    //get touchable->GetCopyNumber(1) for the copy number of the parent volume
    //touchable->GetCopyNumber(2) for the grandparent volume
    //G4VTouchable - a base class for all touchable implementations – defines the following 'requests' (methods) which all touchable have to respond, where depth means always the number of levels up in the tree to be considered:
    //depth = 0 : the bottom level (volume C in B)
    //depth = 1 : the level of its mother volume (volume B in A)
    //depth = 2 : the grandmother volume (volume A in world)
    //G4int G4VTouchable::GetCopyNumber    (    G4int     depth = 0     )     const
        //-> return GetReplicaNumber(depth);
  auto rowNo = touchable->GetReplicaNumber(3);
  auto columnNo = touchable->GetReplicaNumber(2);
  auto layerNumber = touchable->GetReplicaNumber(1);
  G4int hitID = rowNo*fNofPixelsXY*fNofLayers+columnNo*fNofLayers+layerNumber; //iy*fNx*fNz+ix*fNz+iz;
  //auto hitID = rowNo*fNofPixelsXY*(fNofLayers+1)+columnNo*(fNofLayers+1)+layerNumber;
  
  // Get hit accounting data for this cell
  //auto hit = (*fHitsCollection)[layerNumber];
  auto hit = (*fHitsCollection)[hitID];
    
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << hitID;
    G4Exception("CalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }

  // check if it is first touch
  if (hit->GetColumnID()<0) {
    hit->SetRowID(rowNo);
    hit->SetColumnID(columnNo);
    hit->SetLayerID(layerNumber);
    auto depth = touchable->GetHistory()->GetDepth(); //how many level deep in the current tree the volume is
    //G4cout<< " depth "<<depth<< G4endl; // depth = 5
    auto transform = touchable->GetHistory()->GetTransform(depth); //GetTransform(depth) gives the position of the egas or abs cell
    //GetTransform(depth-1) gives the position of the abs+gas cell (bigger cell in the "layer partition")
    transform.Invert();
    hit->SetRot(transform.NetRotation());
    hit->SetPos(transform.NetTranslation());
    hit->SetTime(step->GetPreStepPoint()->GetGlobalTime());
    hit->SetTimeLoc(step->GetPreStepPoint()->GetLocalTime());
    hit->SetParticleID(particlePDG);
  }
    
    
    
  //Get hit for total accounting
  //auto hitIDTotal = rowNo*fNofPixelsXY*(fNofLayers+1)+columnNo*(fNofLayers+1)+fNofLayers; //fNofLayers is the layer ID of the last plane
  //auto hitTotal = (*fHitsCollection)[hitIDTotal];
  //?il totale va tutto in una cella?//
  //auto hitTotal = (*fHitsCollection)[fHitsCollection->entries()-1];

  // Add values
  hit->Add(edep, stepLength);
  //hitTotal->Add(edep, stepLength);
  // add energy deposition
  //hit->AddEdep(edep);


  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) {
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl
       << "-------->Hits Collection: in this event they are " << nofHits
       << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
