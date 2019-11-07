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
// $Id: EventAction.cc 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"
#include <iomanip>

#include <stdio.h>
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* DA, RunAction* RA)
:G4UserEventAction(), fDetectorconstruction(DA), fRunAction(RA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
 //geometry parameters
 Det1SizeZ=fDetectorconstruction->GetSize1();
 Det2SizeZ=fDetectorconstruction->GetSize2();
 Dist=fDetectorconstruction->GetDist();
 Strip1Depth=fDetectorconstruction->GetStrip1Depth();
 Strip1Length=fDetectorconstruction->GetStrip1Length();
 Strip2Depth=fDetectorconstruction->GetStrip2Depth();
 Strip2Length=fDetectorconstruction->GetStrip2Length();
 StripDist=fDetectorconstruction->GetStripDist();
 StripWidth=fDetectorconstruction->GetStripWidth();
 StripPitch = StripDist + StripWidth;
 posEndArm1=-(fDetectorconstruction->Getpos_EndArm1Abs());
 posBeginningArm2=fDetectorconstruction->Getpos_BeginningArm2Abs();
 BPIXSizeZ=fDetectorconstruction->GetSizeBPIX();
 pixelX=fDetectorconstruction->GetPixelPitchX();
 pixelY=fDetectorconstruction->GetPixelPitchY();
 pixelHalf=fDetectorconstruction->GetPixelPitchHalf();
 pixelDepth=fDetectorconstruction->GetPixelDepth();
 XangleDUT=fDetectorconstruction->GetDUTangleX();
 XangleBPIX=fDetectorconstruction->GetBPIXangleX();
 YangleBPIX=fDetectorconstruction->GetBPIXangleY();
 ElField1=fDetectorconstruction->GetElField1();
 ElField2=fDetectorconstruction->GetElField2();
 totalNbStrips = 1016;
 ROChor = 52*pixelX + 2*pixelX;
 ROCvert = 80*pixelY + pixelY;
 fPixThreshold = 1700;

 //G4cout 
     //<< "\n DetectorConstruction parameters assigned. "
     //<< G4endl;

 // initialisation per event
 fEnergyDeposit1  =  fEnergyDeposit2  = fEnergySecondary1 = fEnergySecondary2 = fEnergyTertiary1 = fEnergyTertiary2 = 0.;
 for (G4int j = 0; j <= 1016; j = j + 1) {
    fEnergyStrip1a[j] = 0.;
    fEnergyStrip1b[j] = 0.;
    fEnergyStrip2a[j] = 0.;
    fEnergyStrip2b[j] = 0.;
    fWeightStrip1a[j] = 0.;
    fWeightStrip1b[j] = 0.;
    fWeightStrip2a[j] = 0.;
    fWeightStrip2b[j] = 0.;
    fNbHitsStrip1a[j] = 0;
    fNbHitsStrip1b[j] = 0;
    fNbHitsStrip2a[j] = 0;
    fNbHitsStrip2b[j] = 0;
    fChargeStrip1a[j] = 0;
    fChargeStrip1b[j] = 0;
    fChargeStrip2a[j] = 0;
    fChargeStrip2b[j] = 0;
 }

 totalNbHitStrips = 0;

 for (G4int je1 = 0; je1 <= 2; je1 = je1 + 1) {
    fHitsMultiplicityPix[je1] = 0;
    for (G4int je2 = 0; je2 <= 16; je2 = je2 + 1) {
 	for (G4int je3 = 0; je3 <= 80; je3 = je3 + 1) {
	   for (G4int je4 = 0; je4 <= 52; je4 = je4 + 1) {
		fEnergyPixel[je1][je2][je3][je4] = 0.;
		fChargePixel[je1][je2][je3][je4] = 0;
		fNbHitsPixel[je1][je2][je3][je4] = 0;
 	   }
 	}
    }
 }

 for (G4int je5 = 0; je5 <= 16; je5 = je5 + 1) {
    for (G4int je6 = 0; je6 <= 80; je6 = je6 + 1) {
       for (G4int je7 = 0; je7 <= 52; je7 = je7 + 1) {
          //fWeightPixel1[je5][je6][je7] = 0.;
          fWeightPixel2[je5][je6][je7] = 0.;
       }
    }
 }

 for (G4int je8 = 0; je8 <= 216; je8 = je8 + 1) {
    for (G4int je9 = 0; je9 <= 864; je9 = je9 + 1) {
       fWeightPixel1[je8][je9] = 0.;
       fEnergyPixel1[je8][je9] = 0.;
       fChargePixel1Matrix[je8][je9] = 0.;
    }
 }

 fChargeStrip1 = fChargeStrip2 = 0;
 fChargePixel1 = fChargePixel2 = 0;
 fCS1a = fCS1b = fCS2a = fCS2b = 0;
 fChargePix1 = fChargePix2 = 0;

 fPrimaryTrackLength = fSecondaryTrackLength = fSecondaryDetTrackLength = fTertiaryTrackLength = 0.;
 fSecondaryxPolarization = fSecondaryyPolarization = fSecondaryzPolarization = fTertiaryxPolarization = fTertiaryyPolarization = fTertiaryzPolarization = 0.;
 fHitsMultiplicity1b = fHitsMultiplicity2b = 0;
 fHitSensor1 = fHitSensor2 = 0;
 fHitPixelDet1 = fHitPixelDet2 = 0;

 fMomDir1x = fMomDir1y = fMomDir1z =  fMomDir2x = fMomDir2y = fMomDir2z = 0.;

 fPointPix1ent = fPointPix1entReal = fPointPix1mid = fPointPix1ex = fPointDet1ent = fPointDet1entReal = fPointDet1mid = fPointDet1ex = fPointDet2ent = fPointDet2entReal = fPointDet2mid = fPointDet2ex = fPointPix2ent = fPointPix2entReal = fPointPix2mid = fPointPix2ex = G4ThreeVector(0, 0, 0);

 A1 = A = A2 = At = AReal = B1 = B = B2 = Bt = BReal = C1 = C = C2 = Ct = CReal = D1 = D = D2 = Dt = DReal = E1 = E = E2 = Et = EReal = F1 = F = F2 = Ft = FReal = G1 = G = G2 = Gt = GReal = H1 = H = H2 = Ht = HReal = I1 = I = I2 = It = IReal = J1 = J = J2 = Jt = JReal = Br1 = Br2 = Cr1 = Cr2 = Br1t = Br2t = Cr1t = Cr2t = Br1Real = Br2Real = Cr1Real = Cr2Real = Rav = G4ThreeVector(0, 0, 0);

 fTrack1 = fTrack2 = 0.;

 dB = dC = dBAD = dCAD = dAD = 0.;
 B1Bx = B1By = C1Cx = C1Cy = 0.;

 fDiff1x = fDiff1y = fDiff2x = fDiff2y = 0.;

 xBA = xBD = xCA = xCD = xAD = xBAD = xCAD = B1 = C1 = G4ThreeVector(0, 0, 0);

 //G4cout 
     //<< "\n End of BeginningOfEventAction. "
     //<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 // get event ID
 G4int evtNb = evt->GetEventID();
 //G4cout 
    //<< "\n Event ID = " 
    //<< evtNb
    //<< G4endl;

 fRunAction->AddEnergyDeposit1(fEnergyDeposit1);
 fRunAction->AddEnergyDeposit2(fEnergyDeposit2);

 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 analysisManager->SetVerboseLevel(1);
 analysisManager->FillH1(1, fEnergyDeposit1);
 analysisManager->FillH1(2, fEnergyDeposit2);
 analysisManager->FillH1(3, fEnergySecondary1);
 analysisManager->FillH1(4, fEnergySecondary2);
 analysisManager->FillH1(5, fEnergyTertiary1);
 analysisManager->FillH1(6, fEnergyTertiary2);
 analysisManager->FillH1(7, fEnergyDeposit1+fEnergySecondary1);
 analysisManager->FillH1(8, fEnergyDeposit2+fEnergySecondary2);
 analysisManager->FillH1(21, fSecondaryxPolarization);
 analysisManager->FillH1(22, fSecondaryyPolarization);
 analysisManager->FillH1(23, fSecondaryzPolarization);
 analysisManager->FillH1(24, fTertiaryxPolarization);
 analysisManager->FillH1(25, fTertiaryyPolarization);
 analysisManager->FillH1(26, fTertiaryzPolarization);
 analysisManager->FillH1(27, fPrimaryTrackLength);
 analysisManager->FillH1(28, fSecondaryTrackLength);
 analysisManager->FillH1(29, fTertiaryTrackLength);
 analysisManager->FillH1(30, fSecondaryDetTrackLength);
 analysisManager->FillH1(31, fTrack2);
 for (G4int i = 52; i < 60; i = i + 1) {
   analysisManager->FillH1(i, fEnergyStrip1a[i+453]);
 }
 for (G4int p = 60; p < 68; p = p + 1) {
   analysisManager->FillH1(p, fEnergyStrip1b[p+445]);
 }
 for (G4int l = 68; l < 76; l = l + 1) {
   analysisManager->FillH1(l, fEnergyStrip2a[l+437]);
 }
 for (G4int n = 76; n < 84; n = n + 1) {
   analysisManager->FillH1(n, fEnergyStrip2b[n+429]);
 }

 // fill ntuple and charge per strip histograms
 for (G4int k = 0; k < 1017; k++) {
   analysisManager->FillNtupleDColumn(k, fEnergyStrip1a[k]);
   analysisManager->FillNtupleDColumn(k+1017, fEnergyStrip1b[k]);
   analysisManager->FillNtupleDColumn(k+2*1017, fEnergyStrip2a[k]);
   analysisManager->FillNtupleDColumn(k+3*1017, fEnergyStrip2b[k]);

   //charge collected per strip in number of electrons [e]
   fChargeStrip1a[k] = (fEnergyStrip1a[k])/(3.67*eV);
   fChargeStrip1b[k] = (fEnergyStrip1b[k])/(3.67*eV);
   fChargeStrip2a[k] = (fEnergyStrip2a[k])/(3.67*eV);
   fChargeStrip2b[k] = (fEnergyStrip2b[k])/(3.67*eV);

   fCS1a = fChargeStrip1a[k];
   //G4cout 
     //<< "\n fCS1a = "
     //<< fCS1a;
   fCS1b = fChargeStrip1b[k];
   //G4cout 
     //<< "\n fCS1b = "
     //<< fCS1b;
   fCS2a = fChargeStrip2a[k];
   //G4cout 
     //<< "\n fCS2a = "
     //<< fCS2a;
   fCS2b = fChargeStrip2b[k];
   //G4cout 
     //<< "\n fCS2b = "
     //<< fCS2b;
 
   if (fCS1a >= 5000)  {
      fRunAction->AddNbHitsStrip1a(k);
      fHitSensor1 = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 1 Row: a Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip1a advanced"
        //<< G4endl;

      fChargeStrip1 = fChargeStrip1 + 5000;
      totalNbHitStrips++;
   }
   if (fCS1b >= 5000)  {
      fRunAction->AddNbHitsStrip1b(k);
      fHitSensor1 = 1;
      fHitsMultiplicity1b ++;
      //G4cout 
         //<< "\n Strip hit: Sensor: 1 Row: b Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip1b advanced"
        //<< G4endl;

      fChargeStrip1 = fChargeStrip1 + 5000;
      totalNbHitStrips++;
   }
   if (fCS2a >= 5000)  {
      fRunAction->AddNbHitsStrip2a(k);
      fHitSensor2 = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 2 Row: a Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip2a advanced"
        //<< G4endl;

      fChargeStrip2 = fChargeStrip2 + 5000;
      totalNbHitStrips++;
   }
   if (fCS2b >= 5000)  {
      fRunAction->AddNbHitsStrip2b(k);
      fHitSensor2 = 1;
      fHitsMultiplicity2b ++;
      //G4cout 
         //<< "\n Strip hit: Sensor: 2 Row: b Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip2b advanced"
        //<< G4endl;

      fChargeStrip2 = fChargeStrip2 + 5000;
      totalNbHitStrips++;
   }
 }

 //G4cout 
    //<< "\n totalNbHitStrips = "
    //<< totalNbHitStrips
    //<< G4endl;

 //charge weights and entrance positions of primaries in strip sensors per event
 for (G4int k = 0; k < 1017; k++) {
   fCS1a = fChargeStrip1a[k];
   //G4cout 
     //<< "\n fCS1a = "
     //<< fCS1a;
   fCS1b = fChargeStrip1b[k];
   //G4cout 
     //<< "\n fCS1b = "
     //<< fCS1b;
   fCS2a = fChargeStrip2a[k];
   //G4cout 
     //<< "\n fCS2a = "
     //<< fCS2a;
   fCS2b = fChargeStrip2b[k];
   //G4cout 
     //<< "\n fCS2b = "
     //<< fCS2b;

  if (fCS1a >= 5000)  {
      fWeightStrip1a[k] = 5000/fChargeStrip1;

      fStripCenterNRX = Strip1Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = -Dist/2 - Det1SizeZ;
      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      fPointDet1ent = fPointDet1ent + G4ThreeVector((fStripCenterX*5000/fChargeStrip1), (fStripCenterY*5000/fChargeStrip1), (fStripCenterZ*5000/fChargeStrip1));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS1a = "
         //<< fCS1a;
      //G4cout 
         //<< "\n fChargeStrip1 = "
         //<< fChargeStrip1;
      //G4cout 
         //<< "\n fWeightStrip1a[k] = "
         //<< fWeightStrip1a[k];

      //G4cout 
         //<< "\n fPointDet1ent = "
         //<< fPointDet1ent;
   }
   if (fCS1b >= 5000)  {
      fWeightStrip1b[k] = 5000/fChargeStrip1;

      fStripCenterNRX = -Strip1Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = -Dist/2 - Det1SizeZ;
      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      fPointDet1ent = fPointDet1ent + G4ThreeVector((fStripCenterX*5000/fChargeStrip1), (fStripCenterY*5000/fChargeStrip1), (fStripCenterZ*5000/fChargeStrip1));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS1b = "
         //<< fCS1b;
      //G4cout 
         //<< "\n fChargeStrip1 = "
         //<< fChargeStrip1;
      //G4cout 
         //<< "\n fWeightStrip1b[k] = "
         //<< fWeightStrip1b[k];

      //G4cout 
         //<< "\n fPointDet1ent = "
         //<< fPointDet1ent;
   }
   if (fCS2a >= 5000)  {
      fWeightStrip2a[k] = 5000/fChargeStrip2;

      fStripCenterNRX = Strip2Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = Dist/2;
      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      fPointDet2ent = fPointDet2ent + G4ThreeVector((fStripCenterX*5000/fChargeStrip2), (fStripCenterY*5000/fChargeStrip2), (fStripCenterZ*5000/fChargeStrip2));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS2a = "
         //<< fCS2a;
      //G4cout 
         //<< "\n fChargeStrip2 = "
         //<< fChargeStrip2;
      //G4cout 
         //<< "\n fWeightStrip2a[k] = "
         //<< fWeightStrip2a[k];

      //G4cout 
         //<< "\n fPointDet2ent = "
         //<< fPointDet2ent;
   }
   if (fCS2b >= 5000)  {
      fWeightStrip2b[k] = 5000/fChargeStrip2;

      fStripCenterNRX = -Strip2Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = Dist/2;
      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      fPointDet2ent = fPointDet2ent + G4ThreeVector((fStripCenterX*5000/fChargeStrip2), (fStripCenterY*5000/fChargeStrip2), (fStripCenterZ*5000/fChargeStrip2));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS2b = "
         //<< fCS2b;
      //G4cout 
         //<< "\n fChargeStrip2 = "
         //<< fChargeStrip2;
      //G4cout 
         //<< "\n fWeightStrip2b[k] = "
         //<< fWeightStrip2b[k];

      //G4cout 
         //<< "\n fPointDet2ent = "
         //<< fPointDet2ent;
   }
 }

 if ((fHitSensor1 == 1) && (fHitSensor2 == 1))  {
      fRunAction->AddNbCoinc();
 }

 // fill ntuple and charge per pixel histograms
 G4int kec = 0;
 for (G4int ke1 = 0; ke1 < 3; ke1++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
 	for (G4int ke3 = 0; ke3 < 81; ke3++) {
 	   for (G4int ke4 = 0; ke4 < 53; ke4++) {
		//analysisManager->FillNtupleDColumn(kec+8*1017, fEnergyPixel[ke1][ke2][ke3][ke4]);

   		//charge collected per pixel in number of electrons [e]
   		fChargePixel[ke1][ke2][ke3][ke4] = (fEnergyPixel[ke1][ke2][ke3][ke4])/(3.67*eV);

   		fChargePix = fChargePixel[ke1][ke2][ke3][ke4];
   		//G4cout 
     		  //<< "\n fChargePix = "
     		  //<< fChargePix;
     		  //<< " e" 
                  //<< G4endl;
 
   		if (fChargePix >= fPixThreshold)  {
      		   fRunAction->AddNbHitsPixel(ke1, ke2, ke3, ke4);
                   fHitsMultiplicityPix[ke1] ++;
      		   //G4cout 
        	      //<< "\n NbHitsPixel advanced"
        	      //<< G4endl;
		   kec++;
   		}
 	   }
 	}
    }
 }

 analysisManager->AddNtupleRow();

 analysisManager->FillH1(189, fHitsMultiplicity1b);
 //G4cout 
    //<< "\n hit multiplicity of sensor 1: "
    //<< fHitsMultiplicity1b
    //<< G4endl;
 analysisManager->FillH1(190, fHitsMultiplicity2b);
 //G4cout 
    //<< "\n hit multiplicity of sensor 2: "
    //<< fHitsMultiplicity2b
    //<< G4endl;

 analysisManager->FillH1(197, fHitsMultiplicityPix[1]);
 //G4cout 
    //<< "\n hit multiplicity of BPIX module 1: "
    //<< fHitsMultiplicityPix[1]
    //<< G4endl;
 analysisManager->FillH1(198, fHitsMultiplicityPix[2]);
 //G4cout 
    //<< "\n hit multiplicity of BPIX module 2: "
    //<< fHitsMultiplicityPix[2]
    //<< G4endl;


 fCharge1 = fEnergyDeposit1/(3.67*eV);
 fCharge2 = fEnergyDeposit2/(3.67*eV);
 analysisManager->FillH1(155, fCharge1);
 analysisManager->FillH1(156, fCharge2);

 for (G4int ke2 = 0; ke2 < 17; ke2++) {
   for (G4int ke3 = 0; ke3 < 81; ke3++) {
     for (G4int ke4 = 0; ke4 < 53; ke4++) {
       //charge collected per pixel in number of electrons [e]
       //fChargePixel[1][ke2][ke3][ke4] = (fEnergyPixel[1][ke2][ke3][ke4])/(3.67*eV);
       fChargePixel[2][ke2][ke3][ke4] = (fEnergyPixel[2][ke2][ke3][ke4])/(3.67*eV);

       //fChargePix1 = fChargePixel[1][ke2][ke3][ke4];
       fChargePix2 = fChargePixel[2][ke2][ke3][ke4];
 
       /*if (fChargePix1 >= fPixThreshold)  {
	 fChargePixel1 = fChargePixel1 + fChargePix1;
 	 //G4cout 
    	    //<< "\n Pixel hit: Module: 1 ROC: "
    	    //<< ke2
    	    //<< " row: "
    	    //<< ke3
    	    //<< " col: "
    	    //<< ke4
    	    //<< G4endl;
       }*/
       if (fChargePix2 >= fPixThreshold)  {
	 fChargePixel2 = fChargePixel2 + fChargePix2;
 	 //G4cout 
    	    //<< "\n Pixel hit: Module: 2 ROC: "
    	    //<< ke2
    	    //<< " row: "
    	    //<< ke3
    	    //<< " col: "
    	    //<< ke4
    	    //<< G4endl;
       }
     }
   }
 }

 for (G4int je8 = 0; je8 <= 216; je8 = je8 + 1) {
    for (G4int je9 = 0; je9 <= 864; je9 = je9 + 1) {
       fChargePixel1Matrix[je8][je9] = (fEnergyPixel1[je8][je9])/(3.67*eV);
       fChargePix1 = fChargePixel1Matrix[je8][je9]; 
       if (fChargePix1 >= fPixThreshold)  {
	 fChargePixel1 = fChargePixel1 + fChargePix1;
 	 //G4cout 
    	    //<< "\n fChargePixel 1: "
    	    //<< fChargePixel
    	    //<< " row: "
    	    //<< je8
    	    //<< " col: "
    	    //<< je9
    	    //<< G4endl;
       }
    }
 }


 //charge weights and entrance positions of primaries in BPIX modules per event
 for (G4int je8 = 0; je8 <= 216; je8 = je8 + 1) {
    for (G4int je9 = 0; je9 <= 864; je9 = je9 + 1) {
       fChargePix1 = fChargePixel1Matrix[je8][je9];
       if (fChargePix1 >= fPixThreshold)  {
	 fHitPixelDet1 = 1;
         fWeightPixel1[je8][je9] = fChargePix1/fChargePixel1;
         fPixelCenterNRX = (432 - je9)*pixelHalf + pixelHalf/2;     
     	 fPixelCenterNRY = (108 - je8)*pixelHalf + pixelHalf/2;

	 fPixelCenterNRZ = posEndArm1 - BPIXSizeZ;	 

  	 fPixelCenterNRX1 = fPixelCenterNRX;
  	 fPixelCenterNRY1 = fPixelCenterNRY*cos(XangleBPIX) - fPixelCenterNRZ*sin(XangleBPIX);
   	 fPixelCenterNRZ1 = fPixelCenterNRY*sin(XangleBPIX) + fPixelCenterNRZ*cos(XangleBPIX);

  	 fPixelCenterX = fPixelCenterNRX1*cos(YangleBPIX) + fPixelCenterNRZ1*sin(YangleBPIX);
  	 fPixelCenterY = fPixelCenterNRY1;
  	 fPixelCenterZ = -fPixelCenterNRX1*sin(YangleBPIX) + fPixelCenterNRZ1*cos(YangleBPIX);

         fPointPix1ent = fPointPix1ent + G4ThreeVector((fPixelCenterX*fChargePix1/fChargePixel1), (fPixelCenterY*fChargePix1/fChargePixel1), (fPixelCenterZ*fChargePix1/fChargePixel1));
       }
    }
 }


 for (G4int ke2 = 0; ke2 < 17; ke2++) {
   for (G4int ke3 = 0; ke3 < 81; ke3++) {
     for (G4int ke4 = 0; ke4 < 53; ke4++) {
       //fChargePix1 = fChargePixel[1][ke2][ke3][ke4];
       fChargePix2 = fChargePixel[2][ke2][ke3][ke4];
 
       /* if (fChargePix1 >= fPixThreshold)  {
	 fHitPixelDet1 = 1;
         fWeightPixel1[ke2][ke3][ke4] = fChargePix1/fChargePixel1;

	 if ((ke2 == 1) || (ke2 == 9))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 3*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 3*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 3*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 2) || (ke2 == 10))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 2*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 2*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 2*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 3) || (ke2 == 11))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 1*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 1*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 1*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 4) || (ke2 == 12))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 0*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 0*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 0*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 5) || (ke2 == 13))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -0*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -0*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -0*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 6) || (ke2 == 14))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -1*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -1*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -1*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 7) || (ke2 == 15))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -2*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -2*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -2*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 8) || (ke2 == 16))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -3*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -3*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -3*54*pixelX - 53*pixelX;
	   }
	 }

	 if ((ke2>=1) && (ke2<=8))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = (81 - ke3)*pixelY + pixelY/2;
           }
	 }
	 if ((ke2>=9) && (ke2<=16))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = -80*pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = -(ke3 - 1)*pixelY - pixelY/2;
           }
	 }

	 fPixelCenterNRZ = posEndArm1 - BPIXSizeZ;	 

  	 fPixelCenterNRX1 = fPixelCenterNRX;
  	 fPixelCenterNRY1 = fPixelCenterNRY*cos(XangleBPIX) - fPixelCenterNRZ*sin(XangleBPIX);
   	 fPixelCenterNRZ1 = fPixelCenterNRY*sin(XangleBPIX) + fPixelCenterNRZ*cos(XangleBPIX);

  	 fPixelCenterX = fPixelCenterNRX1*cos(YangleBPIX) + fPixelCenterNRZ1*sin(YangleBPIX);
  	 fPixelCenterY = fPixelCenterNRY1;
  	 fPixelCenterZ = -fPixelCenterNRX1*sin(YangleBPIX) + fPixelCenterNRZ1*cos(YangleBPIX);

         fPointPix1ent = fPointPix1ent + G4ThreeVector((fPixelCenterX*fChargePix1/fChargePixel1), (fPixelCenterY*fChargePix1/fChargePixel1), (fPixelCenterZ*fChargePix1/fChargePixel1));
       }*/

       if (fChargePix2 >= fPixThreshold)  {
	 fHitPixelDet2 = 1;
         fWeightPixel2[ke2][ke3][ke4] = fChargePix2/fChargePixel2;

	 if ((ke2 == 1) || (ke2 == 9))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 3*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 3*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 3*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 2) || (ke2 == 10))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 2*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 2*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 2*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 3) || (ke2 == 11))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 1*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 1*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 1*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 4) || (ke2 == 12))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = 0*54*pixelX + 53*pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = 0*54*pixelX + (53-ke4)*pixelX + pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = 0*54*pixelX + pixelX;
	   }
	 }
	 if ((ke2 == 5) || (ke2 == 13))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -0*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -0*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -0*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 6) || (ke2 == 14))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -1*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -1*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -1*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 7) || (ke2 == 15))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -2*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -2*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -2*54*pixelX - 53*pixelX;
	   }
	 }
	 if ((ke2 == 8) || (ke2 == 16))  {
	   if (ke4 == 1)  {
	     fPixelCenterNRX = -3*54*pixelX - pixelX;
	   }
	   if ((ke4 >= 2) && (ke4<=51))  {
	     fPixelCenterNRX = -3*54*pixelX + -ke4*pixelX - pixelX/2;
	   }
	   if (ke4 == 52)  {
	     fPixelCenterNRX = -3*54*pixelX - 53*pixelX;
	   }
	 }

	 if ((ke2>=1) && (ke2<=8))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = (81 - ke3)*pixelY + pixelY/2;
           }
	 }
	 if ((ke2>=9) && (ke2<=16))  {
	   if (ke3 == 80) {
	     fPixelCenterNRY = -80*pixelY;
           }
	   if ((ke3 >= 1)&&(ke3 <= 79)) {
	     fPixelCenterNRY = -(ke3 - 1)*pixelY - pixelY/2;
           }
	 }

	 fPixelCenterNRZ = posBeginningArm2;	 

  	 fPixelCenterNRX1 = fPixelCenterNRX;
  	 fPixelCenterNRY1 = fPixelCenterNRY*cos(XangleBPIX) - fPixelCenterNRZ*sin(XangleBPIX);
   	 fPixelCenterNRZ1 = fPixelCenterNRY*sin(XangleBPIX) + fPixelCenterNRZ*cos(XangleBPIX);

  	 fPixelCenterX = fPixelCenterNRX1*cos(YangleBPIX) + fPixelCenterNRZ1*sin(YangleBPIX);
  	 fPixelCenterY = fPixelCenterNRY1;
  	 fPixelCenterZ = -fPixelCenterNRX1*sin(YangleBPIX) + fPixelCenterNRZ1*cos(YangleBPIX);

         fPointPix2ent = fPointPix2ent + G4ThreeVector((fPixelCenterX*fChargePix2/fChargePixel2), (fPixelCenterY*fChargePix2/fChargePixel2), (fPixelCenterZ*fChargePix2/fChargePixel2));
       }
     }
   }
 }

 for (G4int i2 = 157; i2 < 165; i2 = i2 + 1) {
   analysisManager->FillH1(i2, fChargeStrip1a[i2+348]);
 }
 for (G4int p2 = 165; p2 < 173; p2 = p2 + 1) {
   analysisManager->FillH1(p2, fChargeStrip1b[p2+340]);
 }
 for (G4int l2 = 173; l2 < 181; l2 = l2 + 1) {
   analysisManager->FillH1(l2, fChargeStrip2a[l2+332]);
 }
 for (G4int n2 = 181; n2 < 189; n2 = n2 + 1) {
   analysisManager->FillH1(n2, fChargeStrip2b[n2+324]);
 }

 ftheta = acos(fMomDir1x*fMomDir2x + fMomDir1y*fMomDir2y + fMomDir1z*fMomDir2z);
 analysisManager->FillH1(116, ftheta);

 //G4cout
    //<< "\n The deflection angle is: "
    //<< ftheta
    //<< " rad"
    //<< G4endl;



 /* //Old residuals

 xBA = fPointDet1ent - fPointPix1ent;
 xBD = fPointDet1ent - fPointPix2ent;
 xCA = fPointDet2ent - fPointPix1ent;
 xCD = fPointDet2ent - fPointPix2ent;
 xAD = fPointPix2ent - fPointPix1ent;

 xBAD = G4ThreeVector((xBA.y())*(xBD.z()) - (xBA.z())*(xBD.y()), (xBA.z())*(xBD.x()) - (xBA.x())*(xBD.z()), (xBA.x())*(xBD.y()) - (xBA.y())*(xBD.x()));
 xCAD = G4ThreeVector((xCA.y())*(xCD.z()) - (xCA.z())*(xCD.y()), (xCA.z())*(xCD.x()) - (xCA.x())*(xCD.z()), (xCA.x())*(xCD.y()) - (xCA.y())*(xCD.x()));

 dBAD = sqrt(sqr(xBAD.x()) + sqr(xBAD.y()) + sqr(xBAD.z()));
 dCAD = sqrt(sqr(xCAD.x()) + sqr(xCAD.y()) + sqr(xCAD.z()));
 dAD = sqrt(sqr(xAD.x()) + sqr(xAD.y()) + sqr(xAD.z()));

 dB = dBAD/dAD;
 dC = dCAD/dAD;

 //G4cout 
    //<< "\n Distance between B and AD: " 
    //<< G4BestUnit(dB, "Length")
    //<< G4endl;
 if (dB == dB)   {
    analysisManager->FillH1(117, dB);
 }

 //G4cout 
    //<< "\n Distance between C and AD: " 
    //<< G4BestUnit(dC, "Length")
    //<< G4endl;
 if (dC == dC)   {
    analysisManager->FillH1(118, dC);
 }

 B1 = G4ThreeVector((xAD.x())*(fPointDet1ent.z() - fPointPix2ent.z())/(xAD.z()) + fPointPix2ent.x(), (xAD.y())*(fPointDet1ent.z() - fPointPix2ent.z())/(xAD.z()) + fPointPix2ent.y(), fPointDet1ent.z());
 B1Bx = fPointDet1ent.x() - B1.x();
 B1By = fPointDet1ent.y() - B1.y();

 if (B1Bx == B1Bx)  {
   analysisManager->FillH1(119, B1Bx);

 //G4cout 
    //<< "\n fPointDet1ent.x(): " 
    //<< G4BestUnit(fPointDet1ent.x(), "Length")
    //<< "\n B1.x(): " 
    //<< G4BestUnit(B1.x(), "Length")
    //<< G4endl;

 }
 if (B1By == B1By)  {
   analysisManager->FillH1(120, B1By);
   if ((fPointDet1ent.y() > -2*StripPitch/3) && (fPointDet1ent.y() < -StripPitch/90))   {
      analysisManager->FillH1(209, B1By);
   }
   if ((fPointDet1ent.y() > -StripPitch/90) && (fPointDet1ent.y() < StripPitch/90))   {
      analysisManager->FillH1(210, B1By);
   }
   if ((fPointDet1ent.y() > StripPitch/90) && (fPointDet1ent.y() < 2*StripPitch/3))   {
      analysisManager->FillH1(211, B1By);
   }
 }

 C1 = G4ThreeVector((xAD.x())*(fPointDet2ent.z() - fPointPix1ent.z())/(xAD.z()) + fPointPix1ent.x(), (xAD.y())*(fPointDet2ent.z() - fPointPix1ent.z())/(xAD.z()) + fPointPix1ent.y(), fPointDet2ent.z());
 C1Cx = fPointDet2ent.x() - C1.x();
 C1Cy = fPointDet2ent.y() - C1.y();
 //G4cout 
    //<< "\n C'Cx: " 
    //<< G4BestUnit(C1Cx, "Length")
    //<< "\n C'Cy: " 
    //<< G4BestUnit(C1Cy, "Length")
    //<< G4endl;

 if (C1Cx == C1Cx)  {
   analysisManager->FillH1(121, C1Cx);
 }
 if (C1Cy == C1Cy)  {
   analysisManager->FillH1(122, C1Cy);
 }


 //G4cout
    //<< "\n A.x = "
    //<< G4BestUnit(fPointPix1ent.x(), "Length")
    //<< "\n A.y = "
    //<< G4BestUnit(fPointPix1ent.y(), "Length")
    //<< "\n A.z = "
    //<< G4BestUnit(fPointPix1ent.z(), "Length")
    //<< "\n B.x = "
    //<< G4BestUnit(fPointDet1ent.x(), "Length")
    //<< "\n B.y = "
    //<< G4BestUnit(fPointDet1ent.y(), "Length")
    //<< "\n B.z = "
    //<< G4BestUnit(fPointDet1ent.z(), "Length")
    //<< "\n B'.x = "
    //<< G4BestUnit(B1.x(), "Length")
    //<< "\n B'.y = "
    //<< G4BestUnit(B1.y(), "Length")
    //<< "\n B'.z = "
    //<< G4BestUnit(B1.z(), "Length")
    //<< "\n C.x = "
    //<< G4BestUnit(fPointDet2ent.x(), "Length")
    //<< "\n C.y = "
    //<< G4BestUnit(fPointDet2ent.y(), "Length")
    //<< "\n C.z = "
    //<< G4BestUnit(fPointDet2ent.z(), "Length")
    //<< "\n C'.x = "
    //<< G4BestUnit(C1.x(), "Length")
    //<< "\n C'.y = "
    //<< G4BestUnit(C1.y(), "Length")
    //<< "\n C'.z = "
    //<< G4BestUnit(C1.z(), "Length")
    //<< "\n D.x = "
    //<< G4BestUnit(fPointPix2ent.x(), "Length")
    //<< "\n D.y = "
    //<< G4BestUnit(fPointPix2ent.y(), "Length")
    //<< "\n D.z = "
    //<< G4BestUnit(fPointPix2ent.z(), "Length")
    //<< "\n \n fPointPix1entReal.x = "
    //<< G4BestUnit(fPointPix1entReal.x(), "Length")
    //<< "\n fPointPix1entReal.y = "
    //<< G4BestUnit(fPointPix1entReal.y(), "Length")
    //<< "\n fPointPix1entReal.z = "
    //<< G4BestUnit(fPointPix1entReal.z(), "Length")
    //<< "\n fPointDet1entReal.x = "
    //<< G4BestUnit(fPointDet1entReal.x(), "Length")
    //<< "\n fPointDet1entReal.y = "
    //<< G4BestUnit(fPointDet1entReal.y(), "Length")
    //<< "\n fPointDet1entReal.z = "
    //<< G4BestUnit(fPointDet1entReal.z(), "Length")
    //<< "\n fPointDet2entReal.x = "
    //<< G4BestUnit(fPointDet2entReal.x(), "Length")
    //<< "\n fPointDet2entReal.y = "
    //<< G4BestUnit(fPointDet2entReal.y(), "Length")
    //<< "\n fPointDet2entReal.z = "
    //<< G4BestUnit(fPointDet2entReal.z(), "Length")
    //<< "\n fPointPix2entReal.x = "
    //<< G4BestUnit(fPointPix2entReal.x(), "Length")
    //<< "\n fPointPix2entReal.y = "
    //<< G4BestUnit(fPointPix2entReal.y(), "Length")
    //<< "\n fPointPix2entReal.z = "
    //<< G4BestUnit(fPointPix2entReal.z(), "Length")
    //<< G4endl; */



 fDiff1x = fPointDet1ent.x() - fPointDet1entReal.x();
 fDiff1y = fPointDet1ent.y() - fPointDet1entReal.y();
 fDiff2x = fPointDet2ent.x() - fPointDet2entReal.x();
 fDiff2y = fPointDet2ent.y() - fPointDet2entReal.y(); 


 //New residuals

 /*Rav = (AReal+DReal)/2;
 At = AReal - Rav;
 Dt = DReal - Rav;*/

 Rav = (fPointPix1ent+fPointPix2ent)/2;
 At = fPointPix1ent - Rav;
 Dt = fPointPix2ent - Rav;

 M[0][0] = ((At.x())*(At.x()) + (Dt.x())*(Dt.x()))/2;
 M[0][1] = M[1][0] = ((At.x())*(At.y()) + (Dt.x())*(Dt.y()))/2;
 M[0][2] = M[2][0] = ((At.x())*(At.z()) + (Dt.x())*(Dt.z()))/2;
 M[1][1] = ((At.y())*(At.y()) + (Dt.y())*(Dt.y()))/2;
 M[1][2] = M[2][1] = ((At.y())*(At.z()) + (Dt.y())*(Dt.z()))/2;
 M[2][2] = ((At.z())*(At.z()) + (Dt.z())*(Dt.z()))/2;

 //G4cout 
    //<< "\n M[0][0]: " 
    //<< M[0][0]
    //<< "\n M[0][1]: " 
    //<< M[0][1]
    //<< "\n M[0][2]: " 
    //<< M[0][2]
    //<< "\n M[1][0]: " 
    //<< M[1][0]
    //<< "\n M[1][1]: " 
    //<< M[1][1]
    //<< "\n M[1][2]: " 
    //<< M[1][2]areal
    //<< "\n M[2][0]: " 
    //<< M[2][0]
    //<< "\n M[2][1]: " 
    //<< M[2][1]
    //<< "\n M[2][2]: " 
    //<< M[2][2]
    //<< G4endl;

 //for (G4int wi = 0; wi < 3; wi = wi + 1) {
    //for (G4int wj = 0; wj < 3; wj = wj + 1) {
        //G4cout 
    	   //<< "\n M[i][j]: " 
           //<< M[wi][wj]
           //<< G4endl;
    //}
 //}


 //Calculation of eigenvalues of symmetric 3X3 matrix M
 p1 = (M[0][1])*(M[0][1]) + (M[0][2])*(M[0][2]) + (M[1][2])*(M[1][2]);
 traceM = M[0][0] + M[1][1] + M[2][2]; //The sum of all diagonal values
 Iu[0][0] = Iu[1][1] = Iu[2][2] = 1;
 Iu[0][1] = Iu[1][0] = Iu[0][2] = Iu[2][0] = Iu[1][2] = Iu[2][1] = 0; //Iu is the identity matrix
 if (p1 == 0) { //M is diagonal
    eig1 = M[0][0];
    eig2 = M[1][1];
    eig3 = M[2][2];
    lambda = eig1;
    if (eig2>lambda)   {
     lambda = eig2;
    }
    if (eig3>lambda)   {
     lambda = eig3;
    }
 }
 else {
    q = traceM/3;
    p2 = (M[0][0] - q)*(M[0][0] - q) + (M[1][1] - q)*(M[1][1] - q) + (M[2][2] - q)*(M[2][2] - q) + 2*p1;
    p = sqrt(p2/6);
    for (G4int bi = 0; bi < 3; bi = bi + 1) {
    	for (G4int bj = 0; bj < 3; bj = bj + 1) {
           if (p>0) {
              Ba[bi][bj] = (M[bi][bj] - q*Iu[bi][bj])/p;
           }
           else {
              Ba[bi][bj] = 0;
           }
    	}
    }
    detBa = (Ba[0][0])*((Ba[1][1])*(Ba[2][2]) - (Ba[1][2])*(Ba[2][1])) - (Ba[0][1])*((Ba[1][0])*(Ba[2][2]) - (Ba[1][2])*(Ba[2][0])) + (Ba[0][2])*((Ba[1][0])*(Ba[2][1]) - (Ba[1][1])*(Ba[2][0]));
    r = detBa/2;

    //In exact arithmetic for a symmetric matrix -1 <= r <= 1 but computation error can leave it slightly outside this range.
    if (r <= -1) {
       phi = 3.14159/3;
    }
    else if (r >= 1) {
       phi = 0;
    }
    else {
       phi = acos(r)/3;
    }

    //The eigenvalues satisfy eig3 <= eig2 <= eig1
    eig1 = q + 2*p*cos(phi);
    eig3 = q + 2*p*cos(phi + 2*3.14159/3);
    eig2 = 3*q - eig1 - eig3; //Since trace(M) = eig1 + eig2 + eig3
    lambda = eig1;
 }

 //Calculation of the eigenvector of maximum eigenvalue lambda with Gaussian elimination
 Ga[0][0] = M[0][0] - lambda;
 Ga[0][1] = Ga[1][0] = M[0][1];
 Ga[0][2] = Ga[2][0] = M[0][2];
 Ga[1][1] = M[1][1] - lambda;
 Ga[1][2] = Ga[2][1] = M[1][2];
 Ga[2][2] = M[2][2] - lambda;

 //G4cout 
    //<< "\n Before the Gaussian elimination: "   	
    //<< "\n Ga[0][0]: " 
    //<< Ga[0][0]
    //<< "\n Ga[0][1]: " 
    //<< Ga[0][1]
    //<< "\n Ga[0][2]: " 
    //<< Ga[0][2]
    //<< "\n Ga[1][0]: " 
    //<< Ga[1][0]
    //<< "\n Ga[1][1]: " 
    //<< Ga[1][1]
    //<< "\n Ga[1][2]: " 
    //<< Ga[1][2]
    //<< "\n Ga[2][0]: " 
    //<< Ga[2][0]
    //<< "\n Ga[2][1]: " 
    //<< Ga[2][1]
    //<< "\n Ga[2][2]: " 
    //<< Ga[2][2]
    //<< G4endl;
 
 G4int h = 0; //Initialization of the pivot row
 G4int k = 0; //Initialization of the pivot column
 G4int i_max;
 G4double maxGa, checkGa, cGa, Ga0, Ga1, Ga2, fa;
 while ((h<=2)&&(k<=2)) {
    //Find the k-th pivot
    i_max = h; //argmax(ii = h ... 2, abs(G[ii,k]))
    maxGa = abs(Ga[h][k]);
    for (G4int ii = h; ii < 3; ii = ii + 1) {
       checkGa = abs(Ga[ii][k]);
       if (checkGa > maxGa) {
	  maxGa = checkGa;
          i_max = ii;
       }
    }

    cGa = Ga[i_max][k];
    if (cGa == 0) { //No pivot in this column, pass to next column
       k = k+1;
    }
    else {
       //Swap rows h, i_max
       Ga0 = Ga[h][0];
       Ga1 = Ga[h][1];
       Ga2 = Ga[h][2];
       Ga[h][0] = Ga[i_max][0];
       Ga[h][1] = Ga[i_max][1];
       Ga[h][2] = Ga[i_max][2];
       Ga[i_max][0] = Ga0;
       Ga[i_max][1] = Ga1;
       Ga[i_max][2] = Ga2;

       //Do for all rows below pivot:
       for (G4int ij = h+1; ij < 3; ij = ij + 1) {
          fa = Ga[ij][k]/Ga[h][k];
          //Fill with zeros the lowest part of pivot column:
          Ga[ij][k] = 0;
          //Do for all remaining elements in current row:
          for (G4int jj = k+1; jj < 3; jj = jj + 1) {
             Ga[ij][jj] = Ga[ij][jj] - fa*Ga[h][jj];
          }
       }

       //Increase pivot row and column
       h = h+1;
       k = k+1;
    }
 }

 //G4cout 
    //<< "\n After the Gaussian elimination: "   	
    //<< "\n Ga[0][0]: " 
    //<< Ga[0][0]
    //<< "\n Ga[0][1]: " 
    //<< Ga[0][1]
    //<< "\n Ga[0][2]: " 
    //<< Ga[0][2]
    //<< "\n Ga[1][0]: " 
    //<< Ga[1][0]
    //<< "\n Ga[1][1]: " 
    //<< Ga[1][1]
    //<< "\n Ga[1][2]: " 
    //<< Ga[1][2]
    //<< "\n Ga[2][0]: " 
    //<< Ga[2][0]
    //<< "\n Ga[2][1]: " 
    //<< Ga[2][1]
    //<< "\n Ga[2][2]: " 
    //<< Ga[2][2]
    //<< G4endl;

 //Matrix3f A;
 //A << M[0][0], M[0][1], M[0][2], M[1][0], M[1][1], M[1][2], M[2][0], M[2][1], M[2][2];

 //SelfAdjointEigenSolver <Matrix3f> es(A);
 //lambda = es.eigenvalues().[0];
 //lambdaVector = es.eigenvectors().col(0);
 //if (es.eigenvalues().[1]>lambda)   {
     //lambda = es.eigenvalues().[1];
     //lambdaVector = es.eigenvectors().col(1);
 //}
 //if (es.eigenvalues().[2]>lambda)   {
     //lambda = es.eigenvalues().[2];
     //lambdaVector = es.eigenvectors().col(2);
 //}

 //lambdaVector definition
 Ga11 = Ga[1][1];
 Ga00 = Ga[0][0];
 lambdaVectorZ = 1;
 if (Ga11==0) {
    lambdaVectorY = 0;
 }
 else {
    lambdaVectorY = -(Ga[1][2])*lambdaVectorZ/Ga11;
 }
 if (Ga00==0) {
    lambdaVectorX = 0;
 }
 else {
    lambdaVectorX = -(lambdaVectorY*Ga[0][1] + lambdaVectorZ*Ga[0][2])/Ga00;
 }
 
 //G4cout 
    //<< "\n lambdaVectorX: " 
    //<< lambdaVectorX
    //<< "\n lambdaVectorY: " 
    //<< lambdaVectorY
    //<< "\n lambdaVectorZ: " 
    //<< lambdaVectorZ
    //<< G4endl;

 //The fitted AJ line is Rav + t*lambdaVector, t real
 tB = (fPointDet1ent.z() - Rav.z())/lambdaVectorZ;
 Bt = G4ThreeVector((Rav.x() + tB*lambdaVectorX), (Rav.y() + tB*lambdaVectorY), fPointDet1ent.z());
 B1Bx = fPointDet1ent.x() - Bt.x();
 B1By = fPointDet1ent.y() - Bt.y();

 if (B1Bx == B1Bx)  {
   analysisManager->FillH1(119, B1Bx);
 }
 if (B1By == B1By)  {
   analysisManager->FillH1(120, B1By);
 }

 tC = (fPointDet2ent.z() - Rav.z())/lambdaVectorZ;
 Ct = G4ThreeVector((Rav.x() + tC*lambdaVectorX), (Rav.y() + tC*lambdaVectorY), fPointDet2ent.z());
 C1Cx = fPointDet2ent.x() - Ct.x();
 C1Cy = fPointDet2ent.y() - Ct.y();
 if (C1Cx == C1Cx)  {
   analysisManager->FillH1(121, C1Cx);
 }
 if (C1Cy == C1Cy)  {
   analysisManager->FillH1(122, C1Cy);
 }





 if (fHitSensor1 == 1)   {
    analysisManager->FillH1(191, fPointDet1ent.x());
    analysisManager->FillH1(192, fPointDet1ent.y());
    analysisManager->FillH1(193, fPointDet1ent.z());
    analysisManager->FillH1(199, fDiff1x);
    analysisManager->FillH1(200, fDiff1y);

    //analysisManager->FillH2(1, fPointDet1ent.y(), fPointDet1entReal.y());
 }
 if (fHitSensor2 == 1)   {
    analysisManager->FillH1(194, fPointDet2ent.x());
    analysisManager->FillH1(195, fPointDet2ent.y());
    analysisManager->FillH1(196, fPointDet2ent.z());
    analysisManager->FillH1(201, fDiff2x);
    analysisManager->FillH1(202, fDiff2y);

    //analysisManager->FillH2(2, fPointDet2ent.y(), fPointDet2entReal.y());
 }
 if (fHitPixelDet1 == 1)   {
    analysisManager->FillH1(203, fPointPix1ent.x());
    analysisManager->FillH1(204, fPointPix1ent.y());
    analysisManager->FillH1(205, fPointPix1ent.z());
 }
 if (fHitPixelDet2 == 1)   {
    analysisManager->FillH1(206, fPointPix2ent.x());
    analysisManager->FillH1(207, fPointPix2ent.y());
    analysisManager->FillH1(208, fPointPix2ent.z());
 }

 //G4cout 
    //<< "\n fPointDet1ent.x() = " 
    //<< G4BestUnit(fPointDet1ent.x(), "Length")
    //<< "\n fPointDet1entReal.x() = " 
    //<< G4BestUnit(fPointDet1entReal.x(), "Length")
    //<< G4endl;

 //G4cout 
     //<< "\n Last track length of secondary created in detector calculated = " 
     //<< G4BestUnit(fTrack2, "Length")
     //<< "\n End of Event. Secondary track length from detector secondaries = " 
     //<< G4BestUnit(fSecondaryDetTrackLength, "Length")
     //<< G4endl;

 //Visualize event if there is a track longer than 1 cm
 //if (fTrack2 > 1.0*cm)  {
     //G4cout 
         //<< "\n fTrack2 = " 
         //<< G4BestUnit(fTrack2, "Length")
         //<< G4endl;
     //G4EventManager* evMan = G4EventManager::GetEventManager();
     //evMan->KeepTheCurrentEvent();
 //}

 //G4cout 
    //<< "\n End of event! "
    //<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

