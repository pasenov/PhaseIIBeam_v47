#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <TMath.h>

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH1F.h"

#define extendedRootTree true

void group()
{
  TFile *f25M = new TFile("proton_100K_25MeV.root");

  TH1F *hMs1 = (TH1F*)f25M->Get("189");
  TH1F *hMs2 = (TH1F*)f25M->Get("190");
  TH1F *hMp1 = (TH1F*)f25M->Get("197");
  TH1F *hMp2 = (TH1F*)f25M->Get("198");

  TFile *f50G = new TFile("proton_100K_50GeV.root");

  TH1F *hGs1 = (TH1F*)f50G->Get("189");
  TH1F *hGs2 = (TH1F*)f50G->Get("190");
  TH1F *hGp1 = (TH1F*)f50G->Get("197");
  TH1F *hGp2 = (TH1F*)f50G->Get("198");
 

/* ******************************************* */
/*  Store Histograms                           */
/* ******************************************* */

 TFile results("group.root","new");

 TCanvas *c1 = new TCanvas("c1", "hit multiplicity for BPIX module 1");
 hMp1->GetXaxis()->SetTitle("number of hits per event");
 hMp1->GetYaxis()->SetTitle(" ");
 hMp1->Write("hit multiplicity for BPIX module 1");
 hMp1->SetLineColor(kBlue);
 hMp1->Draw("HIST B");
 hGp1->GetXaxis()->SetTitle("number of hits per event");
 hGp1->GetYaxis()->SetTitle(" ");
 hGp1->Write("hit multiplicity for BPIX module 1");
 hGp1->SetLineColor(kRed);
 hGp1->Draw("HIST B same");
 TLegend *legend1 = new TLegend(0.1, 0.3, 0.7, 0.9);
 legend1->SetHeader("Legend");
 legend1->AddEntry(hMp1, "25 MeV", "l");
 legend1->AddEntry(hGp1, "50 GeV", "l");
 legend1->Draw();
 c1->Update();
 c1->Write();

 TCanvas *c2 = new TCanvas("c2", "hit multiplicity for BPIX module 2");
 hMp2->GetXaxis()->SetTitle("number of hits per event");
 hMp2->GetYaxis()->SetTitle(" ");
 hMp2->Write("hit multiplicity for BPIX module 2");
 hMp2->SetLineColor(kBlue);
 hMp2->Draw("HIST B");
 hGp2->GetXaxis()->SetTitle("number of hits per event");
 hGp2->GetYaxis()->SetTitle(" ");
 hGp2->Write("hit multiplicity for BPIX module 2");
 hGp2->SetLineColor(kRed);
 hGp2->Draw("HIST B same");
 TLegend *legend2 = new TLegend(0.1, 0.3, 0.7, 0.9);
 legend2->SetHeader("Legend");
 legend2->AddEntry(hMp2, "25 MeV", "l");
 legend2->AddEntry(hGp2, "50 GeV", "l");
 legend2->Draw();
 c2->Update();
 c2->Write();

 TCanvas *c3 = new TCanvas("c3", "hit multiplicity for 2S sensor 1");
 hMs1->GetXaxis()->SetTitle("number of hits per event");
 hMs1->GetYaxis()->SetTitle(" ");
 hMs1->Write("hit multiplicity for 2S sensor 1");
 hMs1->SetLineColor(kBlue);
 hMs1->Draw("HIST B");
 hGs1->GetXaxis()->SetTitle("number of hits per event");
 hGs1->GetYaxis()->SetTitle(" ");
 hGs1->Write("hit multiplicity for 2S sensor 1");
 hGs1->SetLineColor(kRed);
 hGs1->Draw("HIST B same");
 TLegend *legend3 = new TLegend(0.1, 0.3, 0.7, 0.9);
 legend3->SetHeader("Legend");
 legend3->AddEntry(hMs1, "25 MeV", "l");
 legend3->AddEntry(hGs1, "50 GeV", "l");
 legend3->Draw();
 c3->Update();
 c3->Write();

 TCanvas *c4 = new TCanvas("c4", "hit multiplicity for 2S sensor 2");
 hMs2->GetXaxis()->SetTitle("number of hits per event");
 hMs2->GetYaxis()->SetTitle(" ");
 hMs2->Write("hit multiplicity for 2S sensor 2");
 hMs2->SetLineColor(kBlue);
 hMs2->Draw("HIST B");
 hGs2->GetXaxis()->SetTitle("number of hits per event");
 hGs2->GetYaxis()->SetTitle(" ");
 hGs2->Write("hit multiplicity for 2S sensor 2");
 hGs2->SetLineColor(kRed);
 hGs2->Draw("HIST B same");
 TLegend *legend4 = new TLegend(0.1, 0.3, 0.7, 0.9);
 legend4->SetHeader("Legend");
 legend4->AddEntry(hMs2, "25 MeV", "l");
 legend4->AddEntry(hGs2, "50 GeV", "l");
 legend4->Draw();
 c4->Update();
 c4->Write();

 results.Close();
 gROOT->ProcessLine(".q ");
}
