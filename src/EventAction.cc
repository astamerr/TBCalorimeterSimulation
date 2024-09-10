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
/// \file EventAction.cc
/// \brief Implementation of the B4c::EventAction class



#include "EventAction.hh"
#include "CalorimeterSD.hh"
#include "CalorHit.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

namespace B4c
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{
	firstHadInteractionDepth = -999 * CLHEP::m;
	firstHadInteractionTime = 1000 * CLHEP::s;
	leakageEnergy = 0.;
	//HitLayer.push_back(-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorHitsCollection*
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection
    = static_cast<CalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: "
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: "
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
	firstHadInteractionDepth = -999 * CLHEP::m;
	firstHadInteractionTime = 1000 * CLHEP::s;
	leakageEnergy = 0.;
	AbsoHitEdep.clear();
	AbsoHitLayer.clear();
	GapHitEdep.clear();
	GapHitParticle.clear();
	GapHitX.clear();
	GapHitY.clear();
	GapHitZ.clear();
	GapHitLayer.clear();
	GapHitRow.clear();
	GapHitColumn.clear();
	GapHitTime.clear();
	//HitLayer.push_back(-1);
	//HitRow.push_back(-1);
	//HitColumn.push_back(-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // Get hits collections IDs (only once)
  if ( fAbsHCID == -1 ) {
    fAbsHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
    fGapHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  }

  // Get hits collections
  auto absoHC = GetHitsCollection(fAbsHCID, event);
  auto gapHC = GetHitsCollection(fGapHCID, event);

  // Get hit with total values
  auto absoHit = (*absoHC)[absoHC->entries()-1];
  auto gapHit = (*gapHC)[gapHC->entries()-1];

 
  
  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();

  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;

    // DEBUG //
    /*
    G4int nHits = absoHC->entries();
    G4cout << " number of Hits: " << nHits << G4endl;
    for(G4int i=0; i< nHits; i++ ){
      auto absoHit_cell = (*absoHC)[i];
      auto gapHit_cell = (*gapHC)[i];
      //G4cout << " Edep abs layer " << i <<" : "<<  G4BestUnit(absoHit_cell->GetEdep(), "Energy") << G4endl;
      //G4cout << " Edep gap layer " << i <<" : "<<  G4BestUnit(gapHit_cell->GetEdep(), "Energy") << G4endl;
      G4cout << " Event " <<eventID<<" Hit n. " << i
        <<" Row: "<<gapHit_cell->GetRowID()<<" Column: "<<gapHit_cell->GetColumnID()<<" Layer: "<<gapHit_cell->GetLayerID()
        <<" Position: "<<G4BestUnit(gapHit_cell->GetPos(), "Length")
        <<" Global Time: "<<G4BestUnit(gapHit_cell->GetTime(), "Time")
        //<<" Rotation: "<<gapHit_cell->GetRot()
        <<" E. dep. : "<<  G4BestUnit(gapHit_cell->GetEdep(), "Energy") << G4endl;

   } 
   
   */
    PrintEventStatistics(
      absoHit->GetEdep(), absoHit->GetTrackLength(),
      gapHit->GetEdep(), gapHit->GetTrackLength());
}

    // Fill histograms, ntuple
    //
    // get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    
    G4int nHits = absoHC->entries();
    
    G4int NHitsTOT = 0; G4double NHits1 = 0; G4double NHits2 = 0; G4double NHits3 = 0;
    G4double ave_Ar = 30 * eV;
    G4double E_5MIP = 8.31 * keV; //5 volte energia rilasciata da una MIP in 5 mm di Ar (E_MIP = 2MeV*cm2 / g * 1.662 mg/cm3 * 5mm)
    G4double E_15MIP = 24.93 * keV;
    
    //G4cout << " number of Hits: " << nHits << G4endl;
    for(G4int i=0; i< nHits; i++ ){
      auto absoHit_cell = (*absoHC)[i];
      auto gapHit_cell = (*gapHC)[i];

        //filling ntuples with gap info
        if(gapHit_cell->GetLayerID() != -1){
        
        	NHitsTOT++;    
        	//AbsoHitEdep.push_back(absoHit_cell->GetEdep());
        	GapHitEdep.push_back(gapHit_cell->GetEdep());			//ID 3
        	GapHitParticle.push_back(gapHit_cell->GetParticleID());	//ID 4
        	GapHitX.push_back((gapHit_cell->GetPos())[0]); 		//ID 5
        	GapHitY.push_back((gapHit_cell->GetPos())[1]);		//ID 6
        	GapHitZ.push_back((gapHit_cell->GetPos())[2]);		//ID 7
        	GapHitLayer.push_back(gapHit_cell->GetLayerID()); 		//ID 8
        	GapHitRow.push_back(gapHit_cell->GetRowID());			//ID 9
        	GapHitColumn.push_back(gapHit_cell->GetColumnID());		//ID 10
        	GapHitTime.push_back(gapHit_cell->GetTime());			//ID 11
        	
        }
        
        if(gapHit_cell->GetEdep() > ave_Ar && gapHit_cell->GetEdep()< E_5MIP) NHits1++;
        if(gapHit_cell->GetEdep() > E_5MIP && gapHit_cell->GetEdep()< E_15MIP) NHits2++;
        if(gapHit_cell->GetEdep() > E_15MIP) NHits3++;
       /* 
        //filling ntuples with absorber information
       if(absoHit_cell->GetLayerID() != -1) {
		AbsoHitEdep.push_back(absoHit_cell->GetEdep());
		AbsoHitLayer.push_back(absoHit_cell->GetLayerID());
	}
	*/
       /* 
        	AbsHitEdep.push_back(absoHit_cell->GetEdep());			//ID 3
        	AbsHitParticle.push_back(absoHit_cell->GetParticleID());	//ID 4
        	AbsHitX.push_back((absoHit_cell->GetPos())[0]); 		//ID 5
        	AbsHitY.push_back((absoHit_cell->GetPos())[1]);		//ID 6
        	AbsHitZ.push_back((absoHit_cell->GetPos())[2]);		//ID 7
        	AbsHitLayer.push_back(absoHit_cell->GetLayerID()); 		//ID 8
        	AbsHitRow.push_back(absoHit_cell->GetRowID());		//ID 9
        	AbsHitColumn.push_back(absoHit_cell->GetColumnID());		//ID 10
        	AbsHitTime.push_back(absoHit_cell->GetTime());		//ID 11
        
        } */        
        }
        
        //G4cout << "Event " << eventID << " Ntot " << NHitsTOT << " and Egap vector size " << GapHitEdep.size() << G4endl;
        
        analysisManager->FillNtupleIColumn(0, 0, eventID);
        analysisManager->FillNtupleIColumn(0, 1, NHitsTOT);
        analysisManager->FillNtupleIColumn(0, 2, NHits1);
        analysisManager->FillNtupleIColumn(0, 3, NHits2);
        analysisManager->FillNtupleIColumn(0, 4, NHits3);
        analysisManager->FillNtupleDColumn(0, 14, firstHadInteractionDepth);
        analysisManager->FillNtupleDColumn(0, 15, leakageEnergy);
 	
 	//G4cout << " Event " <<eventID<<" Hit n. " << NHitsTOT << " energy leak " << G4BestUnit(leakageEnergy, "Energy") << G4endl;
 
 	 analysisManager->AddNtupleRow(0);	
  // fill histograms
  
  analysisManager->FillH1(0, absoHit->GetEdep());
  analysisManager->FillH1(1, gapHit->GetEdep());
  analysisManager->FillH1(2, absoHit->GetTrackLength());
  analysisManager->FillH1(3, gapHit->GetTrackLength());

  analysisManager->FillH1(4, NHitsTOT);

}
      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
