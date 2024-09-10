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
/// \file Constants.hh
/// \brief Definition of B5 example constants.

#ifndef B5Constants_h
#define B5Constants_h 1

#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace B4c
{

//constexpr G4int kNofHodoscopes1 = 15;
//constexpr G4int kNofHodoscopes2 = 25;
//constexpr G4int kNofChambers = 5;
//constexpr G4int kNofEmColumns = 20;
//constexpr G4int kNofEmRows = 4;
//constexpr G4int kNofEmCells = kNofEmColumns * kNofEmRows;
//constexpr G4int kNofHadColumns = 10;
//constexpr G4int kNofHadRows = 2;
//constexpr G4int kNofHadCells = kNofHadColumns * kNofHadRows;./
constexpr G4int  fNofLayers = 6; //for TB prototype 6-12 -- for Large Prototype 50 
constexpr G4int  fNofPixelsXY = 20; //for TB prototype 20 -- for Large Prototype 100
constexpr G4double cellSizeXY  = 1.*cm; //1.*cm;
constexpr G4double absoThickness = 40.*mm; // for TB prototype 4 cm
constexpr G4double gapThickness =  6.*mm;// for TB prototype 6 mm

constexpr G4double PCBThickness = 1*mm;
constexpr G4double DriftThickness = 1*mm;

//layerThickness = absoThickness + gapThickness;
//calorSizeXY  = cellSizeXY * fNofPixelsXY; //10.*cm; //
//calorThickness = fNofLayers * layerThickness;
//worldSizeXY = 1.2 * calorSizeXY;
//worldSizeZ  = 1.2 * calorThickness;

}

#endif
