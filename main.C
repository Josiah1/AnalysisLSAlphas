//
//   An example tree looper
//
#include "EventReader.h"
#include <string>
#include <iostream>
#include "TChain.h"
#include "TH2D.h"
#include "TFile.h"
#include <fstream>

using namespace std;

EventReader *TR;

int BeginJob();
int EndJob();

int main(int argc, char **argv)
{
    string inputfile;
    string ouputfile;
    string infofile;

    TChain *ReadinChain = new TChain("Event");

    if (argc != 4)
    {
        cout << "Usage:" << endl;
        cout << "    AnalysisLSAlphas inputfile ouputfile infofile" << endl;
        cout << endl;
        return 1;
    }
    else
    {
        inputfile = argv[1];
        ouputfile = argv[2];
        infofile = argv[3];
        ReadinChain->Add(inputfile.c_str());
    }

    // criteria for Bi212-Po212 cascade
    double Bi212_Ep_threshold_low = 0.7;
    double Bi212_Ep_threshold_high = 3;
    double Bi212_Ed_threshold_low = 1;
    double Bi212_Ed_threshold_high = 1.5;
    double Bi212_Tpd_threshold_low = 1e-6;
    double Bi212_Tpd_threshold_high = 3e-6;

    // criteria for Bi214-Po214 cascade
    double Bi214_Ep_threshold_low = 1.5;
    double Bi214_Ep_threshold_high = 3.5;
    double Bi214_Ed_threshold_low = 0.7;
    double Bi214_Ed_threshold_high = 1.2;
    double Bi214_Tpd_threshold_low = 10e-6;
    double Bi214_Tpd_threshold_high = 400e-6;

    // defining the histograms
    TH2D *h2dEpd_Bi212[4];
    TH2D *h2dZR2_Bi212[4];
    TH1D *h1dDistance_Bi212[4];
    TH2D *h2dEpd_Bi214[4];
    TH2D *h2dZR2_Bi214[4];
    TH1D *h1dDistance_Bi214[4];

    for (int i = 0; i < 4; i++)
    {
        h2dEpd_Bi212[i] = new TH2D(Form("h2dEpd_Bi212_%d", i + 1),
                                   Form("h2dEpd_Bi212_%d; Delayed Energy [MeV]; Prompt Energy [MeV]", i + 1),
                                   400, 0, 4, 400, 0, 4);
        h2dZR2_Bi212[i] = new TH2D(Form("h2dZR2_Bi212_%d", i + 1),
                                   Form("h2dZR2_Bi212_%d; Radius^{2} [m^{2}]; Z [m]", i + 1),
                                   400, 0, 25, 400, -5, 5);
        h1dDistance_Bi212[i] = new TH1D(Form("h1dDistance_Bi212_%d", i + 1),
                                        Form("h1dDistance_Bi212_%d; Distance [m]; Events/ 1 mm", i + 1),
                                        5000, 0, 5);

        h2dEpd_Bi214[i] = new TH2D(Form("h2dEpd_Bi214_%d", i + 1),
                                   Form("h2dEpd_Bi214_%d; Delayed Energy [MeV]; Prompt Energy [MeV]", i + 1),
                                   400, 0, 4, 400, 0, 4);
        h2dZR2_Bi214[i] = new TH2D(Form("h2dZR2_Bi214_%d", i + 1),
                                   Form("h2dZR2_Bi214_%d; Radius^{2} [m^{2}]; Z [m]", i + 1),
                                   400, 0, 25, 400, -5, 5);
        h1dDistance_Bi214[i] = new TH1D(Form("h1dDistance_Bi214_%d", i + 1),
                                        Form("h1dDistance_Bi214_%d; Distance [m]; Events/ 1 mm", i + 1),
                                        5000, 0, 5);
    }

    // Load the ReadinChain Tree
    // Get the ReadinChain Tree reader
    TR = new EventReader;
    TR->Init(ReadinChain);

    /* Begin Job */
    /* --------- */
    if (BeginJob() == 0)
    {
        cout << "BeginJob failed" << endl;
        return 0;
    }

    /* The main loop over every stream entries */
    /* --------------------------------------- */
    unsigned int entries = ReadinChain->GetEntries();
    for (unsigned int entry = 0; entry < entries; entry++)
    {
        unsigned int localentry = ReadinChain->LoadTree(entry);
        int ret = TR->GetEntry(localentry);
        if (ret == 0)
        {
            cout << "Warning: Read error" << endl;
        }

        int iad = TR->Det - 1;

        if (TR->Fold == 2)
        {
            if (TR->T2PrevSubEvt[1] > Bi212_Tpd_threshold_low && TR->T2PrevSubEvt[1] < Bi212_Tpd_threshold_high)
            {
                h2dEpd_Bi212[iad]->Fill(TR->E[1], TR->E[0]);
                h2dZR2_Bi212[iad]->Fill(TR->Z[0] / 1e3, (TR->X[0] * TR->X[0] + TR->Y[0] * TR->Y[0]) / 1e6);
                if (TR->E[0] > Bi212_Ep_threshold_low && TR->E[0] < Bi212_Ep_threshold_high && TR->E[1] > Bi212_Ed_threshold_low && TR->E[1] < Bi212_Ed_threshold_high)
                {
                    h1dDistance_Bi212[iad]->Fill(TR->D2First[1]);
                }
            }
            if (TR->T2PrevSubEvt[1] > Bi214_Tpd_threshold_low && TR->T2PrevSubEvt[1] < Bi214_Tpd_threshold_high)
            {
                h2dEpd_Bi214[iad]->Fill(TR->E[1], TR->E[0]);
                h2dZR2_Bi214[iad]->Fill(TR->Z[0] / 1e3, (TR->X[0] * TR->X[0] + TR->Y[0] * TR->Y[0]) / 1e6);
                if (TR->E[0] > Bi214_Ep_threshold_low && TR->E[0] < Bi214_Ep_threshold_high && TR->E[1] > Bi214_Ed_threshold_low && TR->E[1] < Bi214_Ed_threshold_high)
                {
                    h1dDistance_Bi214[iad]->Fill(TR->D2First[1]);
                }
            }
        }
        /***  End of each entry   ***/
    }

    /* End Job */
    /* ------- */
    if (EndJob() == 0)
    {
        cout << "BeginJob failed" << endl;
        return 0;
    }

    TFile *ouFile = new TFile(ouputfile.c_str(), "RECREATE");
    for (int i = 0; i < 4; i++)
    {
        h1dDistance_Bi212[i]->Write();
        h1dDistance_Bi214[i]->Write();
        h2dEpd_Bi212[i]->Write();
        h2dEpd_Bi214[i]->Write();
        h2dZR2_Bi212[i]->Write();
        h2dZR2_Bi214[i]->Write();
    }
    ouFile->Close();

    ofstream outputFile;
    outputFile.open(infofile.c_str());

    if (outputFile.is_open())
    {

        outputFile << "Good Ending." << endl;
        outputFile.close();
    }
    else
    {
        cerr << "错误：无法打开文件进行写入！" << endl;
        return 1;
    }

    return 0;
}

int BeginJob()
{
    return 1;
}

int EndJob()
{
    return 1;
}
