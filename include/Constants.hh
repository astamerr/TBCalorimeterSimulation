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

// Quantities with "1" refer to the first part of the calorimeter (6-8 layers of 20x20 on xy) and those with "2" refer to the second (4 layers of 50x50)

constexpr G4int  fNofLayers1 = 6; //for TB prototype 6-12; 
constexpr G4double absoThickness1 = 20.*mm; //100.*mm;

constexpr G4int  fNofLayers2 = 4; //for TB prototype 6-12; 
constexpr G4double absoThickness2 = 20.*mm; //100.*mm;

constexpr G4int fNofLayers = fNofLayers1 + fNofLayers2;

constexpr G4int  fNofPixelsXY1 = 20; //for TB prototype 20;
constexpr G4int  fNofPixelsXY2 = 50;

constexpr G4double cellSizeXY  = 1.*cm; //1.*cm;

constexpr G4double gapThickness = 6.*mm;//100.*mm;
constexpr G4double PCBThickness = 3*mm;
constexpr G4double DriftThickness = 3*mm;

constexpr G4double airGapThickness = 1*cm;
//layerThickness = absoThickness + gapThickness;
//calorSizeXY  = cellSizeXY * fNofPixelsXY; //10.*cm; //
//calorThickness = fNofLayers * layerThickness;
//worldSizeXY = 1.2 * calorSizeXY;
//worldSizeZ  = 1.2 * calorThickness;

}

#endif
