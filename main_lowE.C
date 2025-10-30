/**
 * @file main_lowE.C
 * @brief Main analysis program for Rn219-Po215 cascade and Po210 single events
 *
 * This program analyzes low-energy nuclear cascade events from Rn219-Po215
 * decay chain and Po210 single alpha events. It processes ROOT tree data,
 * applies selection criteria, and generates histograms for analysis.
 */

#include "EventReader.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TRandom3.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// Global event reader pointer
EventReader *TR;

// Function declarations
int BeginJob();
int EndJob();

int main(int argc, char **argv) {
  // Command line arguments
  string inputfile;
  string ouputfile;
  string infofile;

  // Create TChain for reading ROOT tree
  TChain *ReadinChain = new TChain("Event");

  // Check command line arguments
  if (argc != 4) {
    cout << "Usage:" << endl;
    cout << "    AnalysisLSAlphas inputfile ouputfile infofile" << endl;
    cout << endl;
    return 1;
  } else {
    inputfile = argv[1];
    ouputfile = argv[2];
    infofile = argv[3];
    ReadinChain->Add(inputfile.c_str());
  }

  // Selection criteria for Rn219-Po215 cascade events
  // Rn219 decays to Po215 with alpha emission, Po215 decays with alpha emission
  double Rn219_Ep_threshold_low =
      0.5; // Lower threshold for prompt energy (MeV)
  double Rn219_Ep_threshold_high =
      1.4; // Upper threshold for prompt energy (MeV)
  double Rn219_Ed_threshold_low =
      0.6; // Lower threshold for delayed energy (MeV)
  double Rn219_Ed_threshold_high =
      1.3; // Upper threshold for delayed energy (MeV)
  double Rn219_Tpd_threshold_low =
      1e-3; // Lower threshold for prompt-delayed time (s)
  double Rn219_Tpd_threshold_high =
      4e-3; // Upper threshold for prompt-delayed time (s)

  // Selection criteria for Po210 single alpha events
  // Po210 decays with alpha emission (5.3 MeV alpha)
  double Po210_E_threshold_low = 0.3;  // Lower threshold for Po210 energy (MeV)
  double Po210_E_threshold_high = 0.7; // Upper threshold for Po210 energy (MeV)

  struct EventBrief {
    double E;
    double X, Y, Z;
  };

  map<int, vector<EventBrief>> singlesPool;

  // defining the histograms
  TH2D *h2dEpd_Rn219[4];
  TH2D *h2dEpdBkg_Rn219[4];
  TH2D *h2dEpdSig_Rn219[4];
  TH2D *h2dZR2_Rn219[4];
  TH2D *h2dZR2Bkg_Rn219[4];
  TH2D *h2dZR2Sig_Rn219[4];
  TH1D *h1dDistance_Rn219[4];
  TH1D *h1dDistanceBkg_Rn219[4];
  TH1D *h1dDistanceSig_Rn219[4];

  TH2D *h2dZR2_Po210[4];
  TH1D *h1dE_Po210[4];

  for (int i = 0; i < 4; i++) {
    h2dEpd_Rn219[i] = new TH2D(
        Form("h2dEpd_Rn219_%d", i + 1),
        Form("h2dEpd_Rn219_%d; Delayed Energy [MeV]; Prompt Energy [MeV]",
             i + 1),
        400, 0, 4, 400, 0, 4);
    h2dEpdBkg_Rn219[i] = new TH2D(
        Form("h2dEpdBkg_Rn219_%d", i + 1),
        Form("h2dEpdBkg_Rn219_%d; Delayed Energy [MeV]; Prompt Energy [MeV]",
             i + 1),
        400, 0, 4, 400, 0, 4);
    h2dEpdSig_Rn219[i] = new TH2D(
        Form("h2dEpdSig_Rn219_%d", i + 1),
        Form("h2dEpdSig_Rn219_%d; Delayed Energy [MeV]; Prompt Energy [MeV]",
             i + 1),
        400, 0, 4, 400, 0, 4);
    h2dZR2_Rn219[i] =
        new TH2D(Form("h2dZR2_Rn219_%d", i + 1),
                 Form("h2dZR2_Rn219_%d; Radius^{2} [m^{2}]; Z [m]", i + 1), 400,
                 0, 6, 400, -3, 3);
    h2dZR2Bkg_Rn219[i] =
        new TH2D(Form("h2dZR2Bkg_Rn219_%d", i + 1),
                 Form("h2dZR2Bkg_Rn219_%d; Radius^{2} [m^{2}]; Z [m]", i + 1),
                 400, 0, 6, 400, -3, 3);
    h2dZR2Sig_Rn219[i] =
        new TH2D(Form("h2dZR2Sig_Rn219_%d", i + 1),
                 Form("h2dZR2Sig_Rn219_%d; Radius^{2} [m^{2}]; Z [m]", i + 1),
                 400, 0, 6, 400, -3, 3);
    h1dDistance_Rn219[i] = new TH1D(
        Form("h1dDistance_Rn219_%d", i + 1),
        Form("h1dDistance_Rn219_%d; Distance [m]; Events/ 1 mm", i + 1), 5000,
        0, 5);
    h1dDistanceBkg_Rn219[i] = new TH1D(
        Form("h1dDistanceBkg_Rn219_%d", i + 1),
        Form("h1dDistanceBkg_Rn219_%d; Distance [m]; Events/ 1 mm", i + 1),
        5000, 0, 5);
    h1dDistanceSig_Rn219[i] = new TH1D(
        Form("h1dDistanceSig_Rn219_%d", i + 1),
        Form("h1dDistanceSig_Rn219_%d; Distance [m]; Events/ 1 mm", i + 1),
        5000, 0, 5);

    h2dZR2_Po210[i] =
        new TH2D(Form("h2dZR2_Po210_%d", i + 1),
                 Form("h2dZR2_Po210_%d; Radius^{2} [m^{2}]; Z [m]", i + 1), 400,
                 0, 6, 400, -3, 3);
    h1dE_Po210[i] =
        new TH1D(Form("h1dE_Po210_%d", i + 1),
                 Form("h1dE_Po210_%d; Energy [MeV]; Events/ 0.01 MeV", i + 1),
                 200, 0, 2);
  }

  // Load the ReadinChain Tree
  // Get the ReadinChain Tree reader
  TR = new EventReader;
  TR->Init(ReadinChain);

  /* Begin Job */
  /* --------- */
  if (BeginJob() == 0) {
    cout << "BeginJob failed" << endl;
    return 0;
  }

  for (int iad = 0; iad < 4; iad++) {
    singlesPool[iad].clear();
  }

  // Main event processing loop
  // --------------------------
  unsigned int entries = ReadinChain->GetEntries();
  cout << "Processing " << entries << " events..." << endl;

  for (unsigned int entry = 0; entry < entries; entry++) {
    // Load tree and get entry
    unsigned int localentry = ReadinChain->LoadTree(entry);
    int ret = TR->GetEntry(localentry);
    if (ret == 0) {
      cout << "Warning: Read error at entry " << entry << endl;
      continue;
    }

    // Get detector ID (0-3)
    int iad = TR->Det - 1;

    // Process cascade events (Fold == 2)
    if (TR->Fold == 2) {
      // Rn219-Po215 cascade selection
      if (TR->T2PrevSubEvt[1] > Rn219_Tpd_threshold_low &&
          TR->T2PrevSubEvt[1] < Rn219_Tpd_threshold_high) {
        // Fill energy correlation histogram for events with distance < 500 mm
        if (TR->D2First[1] < 500) {
          h2dEpd_Rn219[iad]->Fill(TR->E[1], TR->E[0]);
        }
        // Fill distance histogram for events within energy thresholds
        if (TR->E[0] > Rn219_Ep_threshold_low &&
            TR->E[0] < Rn219_Ep_threshold_high &&
            TR->E[1] > Rn219_Ed_threshold_low &&
            TR->E[1] < Rn219_Ed_threshold_high) {
          h1dDistance_Rn219[iad]->Fill(TR->D2First[1] / 1e3); // Convert mm to m
        }
        // Fill position histogram for events with all selection criteria
        if (TR->D2First[1] < 500 && TR->E[0] > Rn219_Ep_threshold_low &&
            TR->E[0] < Rn219_Ep_threshold_high &&
            TR->E[1] > Rn219_Ed_threshold_low &&
            TR->E[1] < Rn219_Ed_threshold_high) {
          h2dZR2_Rn219[iad]->Fill((TR->X[0] * TR->X[0] + TR->Y[0] * TR->Y[0]) /
                                      1e6,         // Convert mm² to m²
                                  TR->Z[0] / 1e3); // Convert mm to m
        }
      }
    }
    // Process single events (Fold == 1) for background estimation and Po210
    // analysis
    else if (TR->Fold == 1) {
      // Store single events with time > 1.5 ms for background construction
      if (TR->T2PrevSubEvt[0] > 1.5e-3) {
        EventBrief aEvt;
        aEvt.E = TR->E[0];
        aEvt.X = TR->X[0];
        aEvt.Y = TR->Y[0];
        aEvt.Z = TR->Z[0];
        singlesPool[iad].push_back(aEvt);
      }

      // Fill Po210 energy spectrum for all single events
      h1dE_Po210[iad]->Fill(TR->E[0]);

      // Fill Po210 position histogram for events within energy thresholds
      if (TR->E[0] > Po210_E_threshold_low &&
          TR->E[0] < Po210_E_threshold_high) {
        h2dZR2_Po210[iad]->Fill((TR->X[0] * TR->X[0] + TR->Y[0] * TR->Y[0]) /
                                    1e6,         // Convert mm² to m²
                                TR->Z[0] / 1e3); // Convert mm to m
      }
    }
  }

  for (int iad = 0; iad < 4; iad++) {
    size_t nEvts = singlesPool[iad].size();
    if (nEvts == 0)
      continue;
    for (size_t entry = 0; entry < nEvts / 2 - 1; entry++) {
      int prompt = entry;
      int delayed = entry + nEvts / 2;

      double dx = singlesPool[iad][prompt].X - singlesPool[iad][delayed].X;
      double dy = singlesPool[iad][prompt].Y - singlesPool[iad][delayed].Y;
      double dz = singlesPool[iad][prompt].Z - singlesPool[iad][delayed].Z;
      double dist = sqrt(dx * dx + dy * dy + dz * dz);
      double r2 = singlesPool[iad][prompt].X * singlesPool[iad][prompt].X +
                  singlesPool[iad][prompt].Y * singlesPool[iad][prompt].Y;
      double Ep = singlesPool[iad][prompt].E;
      double Ed = singlesPool[iad][delayed].E;

      double Tpd = gRandom->Uniform(0, 4e-3);

      if (Tpd > Rn219_Tpd_threshold_low && Tpd < Rn219_Tpd_threshold_high) {
        if (dist < 500) {
          h2dEpdBkg_Rn219[iad]->Fill(Ed, Ep);
        }
        if (Ep > Rn219_Ep_threshold_low && Ep < Rn219_Ep_threshold_high &&
            Ed > Rn219_Ed_threshold_low && Ed < Rn219_Ed_threshold_high) {
          h1dDistanceBkg_Rn219[iad]->Fill(dist / 1e3);
        }
        if (dist < 500 && Ep > Rn219_Ep_threshold_low &&
            Ep < Rn219_Ep_threshold_high && Ed > Rn219_Ed_threshold_low &&
            Ed < Rn219_Ed_threshold_high) {
          h2dZR2Bkg_Rn219[iad]->Fill(r2 / 1e6,
                                     singlesPool[iad][delayed].Z / 1e3);
        }
      }
    }
  }

  for (int iad = 0; iad < 4; iad++) {
    size_t nEvts = singlesPool[iad].size();
    if (nEvts == 0)
      continue;
    for (size_t entry = 0; entry < nEvts / 2 - 1; entry++) {
      int delayed = entry;
      int prompt = entry + nEvts / 2;

      double dx = singlesPool[iad][prompt].X - singlesPool[iad][delayed].X;
      double dy = singlesPool[iad][prompt].Y - singlesPool[iad][delayed].Y;
      double dz = singlesPool[iad][prompt].Z - singlesPool[iad][delayed].Z;
      double dist = sqrt(dx * dx + dy * dy + dz * dz);
      double r2 = singlesPool[iad][prompt].X * singlesPool[iad][prompt].X +
                  singlesPool[iad][prompt].Y * singlesPool[iad][prompt].Y;
      double Ep = singlesPool[iad][prompt].E;
      double Ed = singlesPool[iad][delayed].E;

      double Tpd = gRandom->Uniform(0, 4e-3);

      if (Tpd > Rn219_Tpd_threshold_low && Tpd < Rn219_Tpd_threshold_high) {
        if (dist < 500) {
          h2dEpdBkg_Rn219[iad]->Fill(Ed, Ep);
        }
        if (Ep > Rn219_Ep_threshold_low && Ep < Rn219_Ep_threshold_high &&
            Ed > Rn219_Ed_threshold_low && Ed < Rn219_Ed_threshold_high) {
          h1dDistanceBkg_Rn219[iad]->Fill(dist / 1e3);
        }
        if (dist < 500 && Ep > Rn219_Ep_threshold_low &&
            Ep < Rn219_Ep_threshold_high && Ed > Rn219_Ed_threshold_low &&
            Ed < Rn219_Ed_threshold_high) {
          h2dZR2Bkg_Rn219[iad]->Fill(r2 / 1e6,
                                     singlesPool[iad][delayed].Z / 1e3);
        }
      }
    }
  }

  // normalizing background spectra and doing the subtraction
  for (int iad = 0; iad < 4; iad++) {
    double Ncan_212 = h1dDistance_Rn219[iad]->Integral(2000, 5000);
    int Nbkg_212 = h1dDistanceBkg_Rn219[iad]->Integral(2000, 5000);
    if (Nbkg_212 == 0)
      continue;

    double scaleE_212 = Ncan_212 / Nbkg_212;

    h1dDistanceSig_Rn219[iad]->Add(h1dDistance_Rn219[iad],
                                   h1dDistanceBkg_Rn219[iad], 1., -scaleE_212);
    h2dEpdSig_Rn219[iad]->Add(h2dEpd_Rn219[iad], h2dEpdBkg_Rn219[iad], 1.,
                              -scaleE_212);
    h2dZR2Sig_Rn219[iad]->Add(h2dZR2_Rn219[iad], h2dZR2Bkg_Rn219[iad], 1.,
                              -scaleE_212);
  }

  /* End Job */
  /* ------- */
  if (EndJob() == 0) {
    cout << "EndJob failed" << endl;
    return 0;
  }

  TFile *ouFile = new TFile(ouputfile.c_str(), "RECREATE");
  for (int i = 0; i < 4; i++) {
    h1dDistance_Rn219[i]->Write();
    h2dEpd_Rn219[i]->Write();
    h2dZR2_Rn219[i]->Write();

    h1dDistanceBkg_Rn219[i]->Write();
    h2dEpdBkg_Rn219[i]->Write();
    h2dZR2Bkg_Rn219[i]->Write();

    h1dDistanceSig_Rn219[i]->Write();
    h2dEpdSig_Rn219[i]->Write();
    h2dZR2Sig_Rn219[i]->Write();
  }

  for (int i = 0; i < 4; i++) {
    h2dZR2_Po210[i]->Write();
    h1dE_Po210[i]->Write();
  }
  ouFile->Close();

  ofstream outputFile;
  outputFile.open(infofile.c_str());

  if (outputFile.is_open()) {

    outputFile << "Good Ending." << endl;
    outputFile.close();
  } else {
    cerr << "错误：无法打开文件进行写入！" << endl;
    return 1;
  }

  return 0;
}

int BeginJob() { return 1; }

int EndJob() { return 1; }
