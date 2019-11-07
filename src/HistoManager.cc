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
// $Id: HistoManager.cc 72242 2013-07-12 08:44:19Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>
#include <sstream>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  :fFileName("phaseiiv1")
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
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);   //enable inactivation of histograms

  
  // Define histograms start values
  // const G4int kNbofStrips = 1016;
  const G4int kMaxHisto = 212;

  const G4String id[] = { "0", "1", "2", "3" , "4", "5", "6", "7", "8", "9" , "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "100", "101", "102", "103", "104", "105", "106", "107", "108", "109", "110", "111", "112", "113", "114", "115", "116", "117", "118", "119", "120", "121", "122", "123", "124", "125", "126", "127", "128", "129", "130", "131", "132", "133", "134", "135", "136", "137", "138", "139", "140", "141", "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155", "156", "157", "158", "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172", "173", "174", "175", "176", "177", "178", "179", "180", "181", "182", "183", "184", "185", "186", "187", "188", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198", "199", "200", "201", "202", "203", "204", "205", "206", "207", "208", "209", "210", "211" };
  const G4String title[] =
                { "dummy",                                                       //0
                  "continuous energy loss along primary track for Sensor 1 (simulation)",   //1
                  "continuous energy loss along primary track for Sensor 2 (simulation)",   //2
                  "energy from secondaries for Sensor 1 (simulation)",                      //3
                  "energy from secondaries for Sensor 2 (simulation)",                      //4
                  "energy from tertiaries for Sensor 1 (simulation)",                       //5
                  "energy from tertiaries for Sensor 2 (simulation)",                       //6
                  "total energy lost by primaries for Sensor 1 (simulation)",           //7
                  "total energy lost by primaries for Sensor 2 (simulation)",           //8
                  "energy spectrum of secondary e-+ for Sensor 1 (simulation)",             //9
                  "energy spectrum of secondary e-+ for Sensor 2 (simulation)",             //10
                  "energy spectrum of secondary gamma for Sensor 1 (simulation)",           //11
                  "energy spectrum of secondary gamma for Sensor 2 (simulation)",           //12
                  "energy spectrum of tertiary e-+ for Sensor 1 (simulation)",              //13
                  "energy spectrum of tertiary e-+ for Sensor 2 (simulation)",              //14
                  "energy spectrum of tertiary gamma for Sensor 1 (simulation)",            //15
                  "energy spectrum of tertiary gamma for Sensor 2 (simulation)",            //16
                  "step size for primary (simulation)",                                       //17
                  "step size for secondaries (simulation)",                                   //18
                  "step size for tertiaries (simulation)",                                    //19
                  "step size for secondaries created in sensors (simulation)",              //20
                  "x-polarization of secondaries (simulation)",                               //21
                  "y-polarization of secondaries (simulation)",                               //22
                  "z-polarization of secondaries (simulation)",                               //23
                  "x-polarization of tertiaries (simulation)",                                //24
                  "y-polarization of tertiaries (simulation)",                                //25
                  "z-polarization of tertiaries (simulation)",                                //26
                  "track length for primary (simulation)",                                    //27
                  "total track length of secondaries (simulation)",                           //28
                  "total track length of tertiaries (simulation)",                            //29
                  "total track length of secondaries created in sensors (simulation)",      //30
                  "track length of secondaries created in sensors (simulation)",            //31
                  "number of positrons created in Sensor 1 (simulation)",       		 //32
                  "number of positrons created in Sensor 2 (simulation)",       		 //33
                  "number of electrons created in Sensor 1 (simulation)",       		 //34
                  "number of electrons created in Sensor 2 (simulation)",       		 //35
                  "number of pi+ created in Sensor 1 (simulation)",             		 //36
                  "number of pi+ created in Sensor 2 (simulation)",             		 //37
                  "number of pi- created in Sensor 1 (simulation)",              		 //38
                  "number of pi- created in Sensor 2 (simulation)",             		 //39
                  "number of mu+ created in Sensor 1 (simulation)",             		 //40
                  "number of mu+ created in Sensor 2 (simulation)",             		 //41
                  "number of mu- created in Sensor 1 (simulation)",             		 //42
                  "number of mu- created in Sensor 2 (simulation)",              		 //43
		  "number of delta e- created in Sensor 1 that have reached Sensor 2 (simulation)", //44
		  "number of delta e- created in Sensor 2 that have reached Sensor 1 (simulation)", //45
                  "number of e+ created before Sensor 1 (simulation)",               	 //46
                  "number of e- created before Sensor 1 (simulation)",               	 //47
		  "number of e+ created between Sensor 1 and Sensor 2 (simulation)", 	 //48
		  "number of e- created between Sensor 1 and Sensor 2 (simulation)", 	 //49
                  "number of e+ created after Sensor 2 (simulation)",               	 //50
                  "number of e- created after Sensor 2 (simulation)",                	 //51
                  "continuous energy loss along primary track for strip 505 a Sensor 1 (simulation)",   //52
                  "continuous energy loss along primary track for strip 506 a Sensor 1 (simulation)",   //53
                  "continuous energy loss along primary track for strip 507 a Sensor 1 (simulation)",   //54
                  "continuous energy loss along primary track for strip 508 a Sensor 1 (simulation)",   //55
                  "continuous energy loss along primary track for strip 509 a Sensor 1 (simulation)",   //56
                  "continuous energy loss along primary track for strip 510 a Sensor 1 (simulation)",   //57
                  "continuous energy loss along primary track for strip 511 a Sensor 1 (simulation)",   //58
                  "continuous energy loss along primary track for strip 512 a Sensor 1 (simulation)",   //59
                  "continuous energy loss along primary track for strip 505 b Sensor 1 (simulation)",   //60
                  "continuous energy loss along primary track for strip 506 b Sensor 1 (simulation)",   //61
                  "continuous energy loss along primary track for strip 507 b Sensor 1 (simulation)",   //62
                  "continuous energy loss along primary track for strip 508 b Sensor 1 (simulation)",   //63
                  "continuous energy loss along primary track for strip 509 b Sensor 1 (simulation)",   //64
                  "continuous energy loss along primary track for strip 510 b Sensor 1 (simulation)",   //65
                  "continuous energy loss along primary track for strip 511 b Sensor 1 (simulation)",   //66
                  "continuous energy loss along primary track for strip 512 b Sensor 1 (simulation)",   //67
                  "continuous energy loss along primary track for strip 505 a Sensor 2 (simulation)",   //68
                  "continuous energy loss along primary track for strip 506 a Sensor 2 (simulation)",   //69
                  "continuous energy loss along primary track for strip 507 a Sensor 2 (simulation)",   //70
                  "continuous energy loss along primary track for strip 508 a Sensor 2 (simulation)",   //71
                  "continuous energy loss along primary track for strip 509 a Sensor 2 (simulation)",   //72
                  "continuous energy loss along primary track for strip 510 a Sensor 2 (simulation)",   //73
                  "continuous energy loss along primary track for strip 511 a Sensor 2 (simulation)",   //74
                  "continuous energy loss along primary track for strip 512 a Sensor 2 (simulation)",   //75
                  "continuous energy loss along primary track for strip 505 b Sensor 2 (simulation)",   //76
                  "continuous energy loss along primary track for strip 506 b Sensor 2 (simulation)",   //77
                  "continuous energy loss along primary track for strip 507 b Sensor 2 (simulation)",   //78
                  "continuous energy loss along primary track for strip 508 b Sensor 2 (simulation)",   //79
                  "continuous energy loss along primary track for strip 509 b Sensor 2 (simulation)",   //80
                  "continuous energy loss along primary track for strip 510 b Sensor 2 (simulation)",   //81
                  "continuous energy loss along primary track for strip 511 b Sensor 2 (simulation)",   //82
                  "continuous energy loss along primary track for strip 512 b Sensor 2 (simulation)",   //83
                  "number of electrons created in strip 505 a Sensor 1 (simulation)",   //84
                  "number of electrons created in strip 506 a Sensor 1 (simulation)",   //85
                  "number of electrons created in strip 507 a Sensor 1 (simulation)",   //86
                  "number of electrons created in strip 508 a Sensor 1 (simulation)",   //87
                  "number of electrons created in strip 509 a Sensor 1 (simulation)",   //88
                  "number of electrons created in strip 510 a Sensor 1 (simulation)",   //89
                  "number of electrons created in strip 511 a Sensor 1 (simulation)",   //90
                  "number of electrons created in strip 512 a Sensor 1 (simulation)",   //91
                  "number of electrons created in strip 505 b Sensor 1 (simulation)",   //92
                  "number of electrons created in strip 506 b Sensor 1 (simulation)",   //93
                  "number of electrons created in strip 507 b Sensor 1 (simulation)",   //94
                  "number of electrons created in strip 508 b Sensor 1 (simulation)",   //95
                  "number of electrons created in strip 509 b Sensor 1 (simulation)",   //96
                  "number of electrons created in strip 510 b Sensor 1 (simulation)",   //97
                  "number of electrons created in strip 511 b Sensor 1 (simulation)",   //98
                  "number of electrons created in strip 512 b Sensor 1 (simulation)",   //99
                  "number of electrons created in strip 505 a Sensor 2 (simulation)",   //100
                  "number of electrons created in strip 506 a Sensor 2 (simulation)",   //101
                  "number of electrons created in strip 507 a Sensor 2 (simulation)",   //102
                  "number of electrons created in strip 508 a Sensor 2 (simulation)",   //103
                  "number of electrons created in strip 509 a Sensor 2 (simulation)",   //104
                  "number of electrons created in strip 510 a Sensor 2 (simulation)",   //105
                  "number of electrons created in strip 511 a Sensor 2 (simulation)",   //106
                  "number of electrons created in strip 512 a Sensor 2 (simulation)",   //107
                  "number of electrons created in strip 505 b Sensor 2 (simulation)",   //108
                  "number of electrons created in strip 506 b Sensor 2 (simulation)",   //109
                  "number of electrons created in strip 507 b Sensor 2 (simulation)",   //110
                  "number of electrons created in strip 508 b Sensor 2 (simulation)",   //111
                  "number of electrons created in strip 509 b Sensor 2 (simulation)",   //112
                  "number of electrons created in strip 510 b Sensor 2 (simulation)",   //113
                  "number of electrons created in strip 511 b Sensor 2 (simulation)",   //114
                  "number of electrons created in strip 512 b Sensor 2 (simulation)",   //115
                  "angular straggling of primary particle (rad) (simulation)",              //116
                  "distance between midpoint of primary inside Sensor 1 and AD (simulation)",   //117
                  "distance between midpoint of primary inside Sensor 2 and AD (simulation)",   //118
                  "2S first-plane residuals in x-direction (simulation)",   //119
                  "2S first-plane residuals in y-direction (simulation)",   //120
                  "2S second-plane residuals in x-direction (simulation)",   //121
                  "2S second-plane residuals in y-direction (simulation)",   //122
                  "number of hits for strip 505 a Sensor 1 (simulation)",   //123
                  "number of hits for strip 506 a Sensor 1 (simulation)",   //124
                  "number of hits for strip 507 a Sensor 1 (simulation)",   //125
                  "number of hits for strip 508 a Sensor 1 (simulation)",   //126
                  "number of hits for strip 509 a Sensor 1 (simulation)",   //127
                  "number of hits for strip 510 a Sensor 1 (simulation)",   //128
                  "number of hits for strip 511 a Sensor 1 (simulation)",   //129
                  "number of hits for strip 512 a Sensor 1 (simulation)",   //130
                  "number of hits for strip 505 b Sensor 1 (simulation)",   //131
                  "number of hits for strip 506 b Sensor 1 (simulation)",   //132
                  "number of hits for strip 507 b Sensor 1 (simulation)",   //133
                  "number of hits for strip 508 b Sensor 1 (simulation)",   //134
                  "number of hits for strip 509 b Sensor 1 (simulation)",   //135
                  "number of hits for strip 510 b Sensor 1 (simulation)",   //136
                  "number of hits for strip 511 b Sensor 1 (simulation)",   //137
                  "number of hits for strip 512 b Sensor 1 (simulation)",   //138
                  "number of hits for strip 505 a Sensor 2 (simulation)",   //139
                  "number of hits for strip 506 a Sensor 2 (simulation)",   //140
                  "number of hits for strip 507 a Sensor 2 (simulation)",   //141
                  "number of hits for strip 508 a Sensor 2 (simulation)",   //142
                  "number of hits for strip 509 a Sensor 2 (simulation)",   //143
                  "number of hits for strip 510 a Sensor 2 (simulation)",   //144
                  "number of hits for strip 511 a Sensor 2 (simulation)",   //145
                  "number of hits for strip 512 a Sensor 2 (simulation)",   //146
                  "number of hits for strip 505 b Sensor 2 (simulation)",   //147
                  "number of hits for strip 506 b Sensor 2 (simulation)",   //148
                  "number of hits for strip 507 b Sensor 2 (simulation)",   //149
                  "number of hits for strip 508 b Sensor 2 (simulation)",   //150
                  "number of hits for strip 509 b Sensor 2 (simulation)",   //151
                  "number of hits for strip 510 b Sensor 2 (simulation)",   //152
                  "number of hits for strip 511 b Sensor 2 (simulation)",   //153
                  "number of hits for strip 512 b Sensor 2 (simulation)",   //154
                  "charge in Sensor 1 corresponding to the deposited energy in Sensor 1 (simulation)",   //155
                  "charge in Sensor 2 corresponding to the deposited energy in Sensor 2 (simulation)",   //156
                  "charge in strip 505 a Sensor 1 (simulation)",   //157
                  "charge in strip 506 a Sensor 1 (simulation)",   //158
                  "charge in strip 507 a Sensor 1 (simulation)",   //159
                  "charge in strip 508 a Sensor 1 (simulation)",   //160
                  "charge in strip 509 a Sensor 1 (simulation)",   //161
                  "charge in strip 510 a Sensor 1 (simulation)",   //162
                  "charge in strip 511 a Sensor 1 (simulation)",   //163
                  "charge in strip 512 a Sensor 1 (simulation)",   //164
                  "charge in strip 505 b Sensor 1 (simulation)",   //165
                  "charge in strip 506 b Sensor 1 (simulation)",   //166
                  "charge in strip 507 b Sensor 1 (simulation)",   //167
                  "charge in strip 508 b Sensor 1 (simulation)",   //168
                  "charge in strip 509 b Sensor 1 (simulation)",   //169
                  "charge in strip 510 b Sensor 1 (simulation)",   //170
                  "charge in strip 511 b Sensor 1 (simulation)",   //171
                  "charge in strip 512 b Sensor 1 (simulation)",   //172
                  "charge in strip 505 a Sensor 2 (simulation)",   //173
                  "charge in strip 506 a Sensor 2 (simulation)",   //174
                  "charge in strip 507 a Sensor 2 (simulation)",   //175
                  "charge in strip 508 a Sensor 2 (simulation)",   //176
                  "charge in strip 509 a Sensor 2 (simulation)",   //177
                  "charge in strip 510 a Sensor 2 (simulation)",   //178
                  "charge in strip 511 a Sensor 2 (simulation)",   //179
                  "charge in strip 512 a Sensor 2 (simulation)",   //180
                  "charge in strip 505 b Sensor 2 (simulation)",   //181
                  "charge in strip 506 b Sensor 2 (simulation)",   //182
                  "charge in strip 507 b Sensor 2 (simulation)",   //183
                  "charge in strip 508 b Sensor 2 (simulation)",   //184
                  "charge in strip 509 b Sensor 2 (simulation)",   //185
                  "charge in strip 510 b Sensor 2 (simulation)",   //186
                  "charge in strip 511 b Sensor 2 (simulation)",   //187
                  "charge in strip 512 b Sensor 2 (simulation)",   //188
                  "cluster size per event: Sensor 1 (simulation)",   //189
                  "cluster size per event: Sensor 2 (simulation)",   //190
                  "reconstructed x-position of primary in sensor 1 (simulation)",    //191
                  "reconstructed y-position of primary in sensor 1 (simulation)",    //192
                  "reconstructed z-position of primary in sensor 1 (simulation)",    //193
                  "reconstructed x-position of primary in sensor 2 (simulation)",    //194
                  "reconstructed y-position of primary in sensor 2 (simulation)",    //195
                  "reconstructed z-position of primary in sensor 2 (simulation)",    //196
                  "cluster size per event: BPIX module 1 (simulation)",    //197
                  "cluster size per event: BPIX module 2 (simulation)",    //198
                  "x-value of difference: reconstructed - extrapolated points for the 2S 1st plane (simulation)",     //199
                  "y-value of difference: reconstructed - extrapolated points for the 2S 1st plane (simulation)",     //200
                  "x-value of difference: reconstructed - extrapolated points for the 2S 2nd plane (simulation)",     //201
                  "y-value of difference: reconstructed - extrapolated points for the 2S 2nd plane (simulation)",     //202
                  "reconstructed x-position of primary in pixel module 1 (simulation)",     //203
                  "reconstructed y-position of primary in pixel module 1 (simulation)",     //204
                  "reconstructed z-position of primary in pixel module 1 (simulation)",     //205
                  "reconstructed x-position of primary in pixel module 2 (simulation)",     //206
                  "reconstructed y-position of primary in pixel module 2 (simulation)",     //207
                  "reconstructed z-position of primary in pixel module 2 (simulation)",     //208
                  //"midpoint of primary particle inside Pixel Detector 1 (simulation)",
                  //"midpoint of primary particle inside Strip Detector 1 (simulation)",
                  //"midpoint of primary particle inside Strip Detector 2 (simulation)",
                  //"midpoint of primary particle inside Pixel Detector 2"
		  "2S first-plane residuals in y-direction when reconstructed y = -45 um, single hit in strip 508 (simulation)",     //209
		  "2S first-plane residuals in y-direction when reconstructed y = 0 um, double hit in stips 508, 509 (simulation)",     //210
		  "2S first-plane residuals in y-direction when reconstructed y = 45 um, single hit in strip 509 (simulation)"      //211
                 };
            
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = -0.5;
  G4double vmax = 100.5;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }

  G4int ih0 = analysisManager->CreateH2("0", "dummy", nbins, vmin, vmax,  //real point
						      nbins, vmin, vmax);   //reconstructed point
  analysisManager->SetH2Activation(ih0, false);

  //G4int ih2a = analysisManager->CreateH2("1", "2D histogram of reconstructed vs. real impact point y-coordinate in Sensor 1", nbins, vmin, vmax,  //real point
																		//nbins, vmin, vmax);   //reconstructed point
  //analysisManager->SetH2Activation(ih2a, false);
  //G4int ih2b = analysisManager->CreateH2("2", "2D histogram of reconstructed vs. real impact point y-coordinate in Sensor 2", nbins, vmin, vmax,  //real point
																		//nbins, vmin, vmax);   //reconstructed point
  //analysisManager->SetH2Activation(ih2b, false);

  //Create H3 histograms
  //G4int ih1 = analysisManager->CreateH3("85", "midpoint of primary particle inside Pixel Detector 1", 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, "cm", "cm", "cm");
  //analysisManager->SetH3Activation(ih1, false);
  //G4int ih2 = analysisManager->CreateH3("86", "midpoint of primary particle inside Strip Detector 1", 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, "mm", "mm", "mm");
  //analysisManager->SetH3Activation(ih2, false);
  //G4int ih3 = analysisManager->CreateH3("87", "midpoint of primary particle inside Strip Detector 2", 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, 1000, -2*mm, 2*mm, "mm", "mm", "mm");
  //analysisManager->SetH3Activation(ih3, false);
  //G4int ih4 = analysisManager->CreateH3("88", "midpoint of primary particle inside Pixel Detector 2", 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, 10000, -2*cm, 2*cm, "cm", "cm", "cm");
  //analysisManager->SetH3Activation(ih4, false);

  // nTuples, column l: 1 to 1016 correspond to Detector 1, 2 to 2032 correspond to Detector 2
  //
  analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetFirstNtupleId(1);       
  analysisManager->CreateNtuple("1", "Primary Particle Tuple 1");
  for (G4int l=0; l<1017; l++) {
     std::ostringstream os1;
     os1 <<"EdepStripSensor1a_" << l;
     analysisManager->CreateNtupleDColumn(os1.str());          //
  }
  for (G4int m1=0; m1<1017; m1++) {
     std::ostringstream os2;
     os2 <<"EdepStripSensor1b_" << m1;
     analysisManager->CreateNtupleDColumn(os2.str());          //
  }
  for (G4int ma=0; ma<1017; ma++) {
     std::ostringstream os3;
     os3 <<"EdepStripSensor2a_" << ma;
     analysisManager->CreateNtupleDColumn(os3.str());          //
  }
  for (G4int mb=0; mb<1017; mb++) {
     std::ostringstream os4;
     os4 <<"EdepStripSensor2b_" << mb;
     analysisManager->CreateNtupleDColumn(os4.str());          //
  }

  for (G4int l1=0; l1<1017; l1++) {
     std::ostringstream os5;
     os5 <<"NbHitsSensor1a_" << l1;
     analysisManager->CreateNtupleDColumn(os5.str());          //
  }
  for (G4int l2=0; l2<1017; l2++) {
     std::ostringstream os6;
     os6 <<"NbHitsSensor1b_" << l2;
     analysisManager->CreateNtupleDColumn(os6.str());          //
  }
  for (G4int l3=0; l3<1017; l3++) {
     std::ostringstream os7;
     os7 <<"NbHitsSensor2a_" << l3;
     analysisManager->CreateNtupleDColumn(os7.str());          //
  }
  for (G4int l4=0; l4<1017; l4++) {
     std::ostringstream os8;
     os8 <<"NbHitsSensor2b_" << l4;
     analysisManager->CreateNtupleDColumn(os8.str());          //
  }

  G4int NbBPIX = 2;
  G4int NbROC = 16;
  G4int NbRows = 80;
  G4int NbCols = 52;

  for (G4int ne1=0; ne1<=NbBPIX; ne1++) {
     for (G4int ne2=0; ne2<=NbROC; ne2++) {
        for (G4int ne3=0; ne3<=NbRows; ne3++) {
           for (G4int ne4=0; ne4<=NbCols; ne4++) {
              std::ostringstream os9;
     	      os9 <<"EdepPixel_" << ne1 <<"_" << ne2 <<"_" << ne3 << "_" << ne4;
     	      //analysisManager->CreateNtupleDColumn(os9.str());  //Name convention: EdepPixel_NbBPIX_NbROC_NbRow_NbCol, zeros are dummy
           }
        }
     }
  }

  for (G4int n1=0; n1<=NbBPIX; n1++) {
     for (G4int n2=0; n2<=NbROC; n2++) {
        for (G4int n3=0; n3<=NbRows; n3++) {
           for (G4int n4=0; n4<=NbCols; n4++) {
              std::ostringstream os10;
     	      os10 <<"NbHitsPixel_" << n1 <<"_" << n2 <<"_" << n3 << "_" << n4;
     	      //analysisManager->CreateNtupleDColumn(os10.str());  //Name convention: NbHitsPixel_NbBPIX_NbROC_NbRow_NbCol, zeros are dummy
           }
        }
     }
  }

  //for (G4int n=0; n<1017; n++) {
     //std::ostringstream os3;
     //os3 <<"NbElectronsStripSensor1_" << n;
     //analysisManager->CreateNtupleDColumn(os3.str());          //
  //}
  //for (G4int o=0; o<1017; o++) {
     //std::ostringstream os4;
     //os4 <<"NbElectronsStripSensor2_" << o;
     //analysisManager->CreateNtupleDColumn(os4.str());          //
  //}
  analysisManager->FinishNtuple();
  
  //analysisManager->SetNtupleActivation(false);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
