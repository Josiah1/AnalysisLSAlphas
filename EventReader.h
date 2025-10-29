//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 29 17:55:47 2025 by ROOT version 6.36.04
// from TTree Event/Analysis Tree
// found on file: /afs/ihep.ac.cn/users/l/lijj16/WORK/DATA/Retw/IBD/H_IBD_0p7MeV/Output/EH1/21221.TWin.root
//////////////////////////////////////////////////////////

#ifndef EventReader_h
#define EventReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

class EventReader : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Det;
   Double_t        T2PrevMu;
   Double_t        T2PrevPlMu;
   Double_t        T2PrevIwsMu;
   Double_t        T2PrevAdMu;
   Double_t        T2PrevShMu;
   Int_t           PrevIwsNhits;
   Double_t        T2PrevNetMu;
   Int_t           Fold;
   Int_t           TrigNum[5];   //[Fold]
   Int_t           TrigSec[5];   //[Fold]
   Int_t           TrigNano[5];   //[Fold]
   Double_t        E[5];   //[Fold]
   Double_t        X[5];   //[Fold]
   Double_t        Y[5];   //[Fold]
   Double_t        Z[5];   //[Fold]
   Double_t        D2First[5];   //[Fold]
   Double_t        T2PrevSubEvt[5];   //[Fold]
   Double_t        Quadrant[5];   //[Fold]
   Double_t        Q1[5];   //[Fold]
   Double_t        Q2[5];   //[Fold]
   Double_t        QMax2Sum[5];   //[Fold]

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Det;   //!
   TBranch        *b_T2PrevMu;   //!
   TBranch        *b_T2PrevPlMu;   //!
   TBranch        *b_T2PrevIwsMu;   //!
   TBranch        *b_T2PrevAdMu;   //!
   TBranch        *b_T2PrevShMu;   //!
   TBranch        *b_PrevIwsNhits;   //!
   TBranch        *b_T2PrevNetMu;   //!
   TBranch        *b_Fold;   //!
   TBranch        *b_TrigNum;   //!
   TBranch        *b_TrigSec;   //!
   TBranch        *b_TrigNano;   //!
   TBranch        *b_E;   //!
   TBranch        *b_X;   //!
   TBranch        *b_Y;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_D2First;   //!
   TBranch        *b_T2PrevSubEvt;   //!
   TBranch        *b_Quadrant;   //!
   TBranch        *b_Q1;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_QMax2Sum;   //!

   EventReader(TTree * /*tree*/ =0) : fChain(0) { }
   ~EventReader() override { }
   Int_t  Version() const override { return 2; }
   void   Begin(TTree *tree) override;
   void   SlaveBegin(TTree *tree) override;
   void   Init(TTree *tree) override;
   bool   Notify() override;
   bool   Process(Long64_t entry) override;
   Int_t  GetEntry(Long64_t entry, Int_t getall = 0) override { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   void   SetOption(const char *option) override { fOption = option; }
   void   SetObject(TObject *obj) override { fObject = obj; }
   void   SetInputList(TList *input) override { fInput = input; }
   TList* GetOutputList() const override { return fOutput; }
   void   SlaveTerminate() override;
   void   Terminate() override;

   ClassDefOverride(EventReader,0);
};

#endif

#ifdef EventReader_cxx
void EventReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Det", &Det, &b_Det);
   fChain->SetBranchAddress("T2PrevMu", &T2PrevMu, &b_T2PrevMu);
   fChain->SetBranchAddress("T2PrevPlMu", &T2PrevPlMu, &b_T2PrevPlMu);
   fChain->SetBranchAddress("T2PrevIwsMu", &T2PrevIwsMu, &b_T2PrevIwsMu);
   fChain->SetBranchAddress("T2PrevAdMu", &T2PrevAdMu, &b_T2PrevAdMu);
   fChain->SetBranchAddress("T2PrevShMu", &T2PrevShMu, &b_T2PrevShMu);
   fChain->SetBranchAddress("PrevIwsNhits", &PrevIwsNhits, &b_PrevIwsNhits);
   fChain->SetBranchAddress("T2PrevNetMu", &T2PrevNetMu, &b_T2PrevNetMu);
   fChain->SetBranchAddress("Fold", &Fold, &b_Fold);
   fChain->SetBranchAddress("TrigNum", TrigNum, &b_TrigNum);
   fChain->SetBranchAddress("TrigSec", TrigSec, &b_TrigSec);
   fChain->SetBranchAddress("TrigNano", TrigNano, &b_TrigNano);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("X", X, &b_X);
   fChain->SetBranchAddress("Y", Y, &b_Y);
   fChain->SetBranchAddress("Z", Z, &b_Z);
   fChain->SetBranchAddress("D2First", D2First, &b_D2First);
   fChain->SetBranchAddress("T2PrevSubEvt", T2PrevSubEvt, &b_T2PrevSubEvt);
   fChain->SetBranchAddress("Quadrant", Quadrant, &b_Quadrant);
   fChain->SetBranchAddress("Q1", Q1, &b_Q1);
   fChain->SetBranchAddress("Q2", Q2, &b_Q2);
   fChain->SetBranchAddress("QMax2Sum", QMax2Sum, &b_QMax2Sum);
}

bool EventReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

#endif // #ifdef EventReader_cxx
