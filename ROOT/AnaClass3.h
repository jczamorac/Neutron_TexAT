//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 28 18:29:24 2022 by ROOT version 6.24/06
// from TChain SiTelescope/
//////////////////////////////////////////////////////////

#ifndef AnaClass3_h
#define AnaClass3_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include "TF1.h"
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TH2Poly.h>
#include "TEnv.h"

#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <Math/Vector3D.h>
#include <TMath.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <numeric> // std::iota
// Header file for the classes stored in the TTree if any.

class AnaClass3 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         E0;
   Float_t         E1;
   Float_t         E2;
   Float_t         E3;
   Float_t         E4;
   Float_t         E5;
   Float_t         E6;
   Float_t         E7;
   Float_t         EKin0;
   Float_t         EKin1;
   Float_t         EKin2;
   Float_t         EKin3;
   Float_t         EKin4;
   Float_t         EKin5;
   Float_t         EKin6;
   Float_t         EKin7;
   Float_t         pos_x_det0;
   Float_t         pos_y_det0;
   Float_t         pos_z_det0;
   Float_t         pos_x_det1;
   Float_t         pos_y_det1;
   Float_t         pos_z_det1;
   Float_t         pos_x_det2;
   Float_t         pos_y_det2;
   Float_t         pos_z_det2;
   Float_t         pos_x_det3;
   Float_t         pos_y_det3;
   Float_t         pos_z_det3;
   Float_t         pos_x_det4;
   Float_t         pos_y_det4;
   Float_t         pos_z_det4;
   Float_t         pos_x_det5;
   Float_t         pos_y_det5;
   Float_t         pos_z_det5;
   Float_t         pos_x_det6;
   Float_t         pos_y_det6;
   Float_t         pos_z_det6;
   Float_t         pos_x_det7;
   Float_t         pos_y_det7;
   Float_t         pos_z_det7;
   Float_t         pos_x_target;
   Float_t         pos_y_target;
   Float_t         pos_z_target;
   Float_t         pos_x_dummy;
   Float_t         pos_y_dummy;
   Float_t         pos_z_dummy;
   Float_t         eloss_dummy;
   Float_t         inang_target;
   Float_t         Reaang_target;
   Float_t         ener_det0;
   Float_t         ener_det1;
   Float_t         ener_det2;
   Float_t         ener_det3;
   Float_t         ener_det4;
   Float_t         ener_det5;
   Float_t         ener_det6;
   Float_t         ener_det7;
   Float_t         rThetaCM0;
   Float_t         rThetaCM1;
   Float_t         rThetaCM2;
   Float_t         rThetaCM3;
   Float_t         rThetaCM4;
   Float_t         rThetaCM5;
   Float_t         rThetaCM6;
   Float_t         rThetaCM7;
   Float_t         t_sili0;
   Float_t         t_sili1;
   Float_t         t_sili2;
   Float_t         t_sili3;
   Float_t         t_sili4;
   Float_t         t_sili5;
   Float_t         t_sili6;
   Float_t         t_sili7;
   Float_t         px_dssd;
   Float_t         py_dssd;
   Float_t         pz_dssd;
   Int_t           Strip_Number0;
   Int_t           Strip_Number1;
   Int_t           Strip_Number2;
   Int_t           Strip_Number3;
   Int_t           Strip_Number4;
   Int_t           Strip_Number5;
   Int_t           Strip_Number6;
   Int_t           Strip_Number7;
   Float_t         truthPosx;
   Float_t         truthPosy;
   Float_t         truthPosz;
   Float_t         truthAngle_theta;
   Float_t         truthAngle_phi;
   Float_t         t_dssd;
   Float_t         t_dssd2;
   Float_t         Etot;

   // List of branches
   TBranch        *b_E0;   //!
   TBranch        *b_E1;   //!
   TBranch        *b_E2;   //!
   TBranch        *b_E3;   //!
   TBranch        *b_E4;   //!
   TBranch        *b_E5;   //!
   TBranch        *b_E6;   //!
   TBranch        *b_E7;   //!
   TBranch        *b_EKin0;   //!
   TBranch        *b_EKin1;   //!
   TBranch        *b_EKin2;   //!
   TBranch        *b_EKin3;   //!
   TBranch        *b_EKin4;   //!
   TBranch        *b_EKin5;   //!
   TBranch        *b_EKin6;   //!
   TBranch        *b_EKin7;   //!
   TBranch        *b_pos_x_det0;   //!
   TBranch        *b_pos_y_det0;   //!
   TBranch        *b_pos_z_det0;   //!
   TBranch        *b_pos_x_det1;   //!
   TBranch        *b_pos_y_det1;   //!
   TBranch        *b_pos_z_det1;   //!
   TBranch        *b_pos_x_det2;   //!
   TBranch        *b_pos_y_det2;   //!
   TBranch        *b_pos_z_det2;   //!
   TBranch        *b_pos_x_det3;   //!
   TBranch        *b_pos_y_det3;   //!
   TBranch        *b_pos_z_det3;   //!
   TBranch        *b_pos_x_det4;   //!
   TBranch        *b_pos_y_det4;   //!
   TBranch        *b_pos_z_det4;   //!
   TBranch        *b_pos_x_det5;   //!
   TBranch        *b_pos_y_det5;   //!
   TBranch        *b_pos_z_det5;   //!
   TBranch        *b_pos_x_det6;   //!
   TBranch        *b_pos_y_det6;   //!
   TBranch        *b_pos_z_det6;   //!
   TBranch        *b_pos_x_det7;   //!
   TBranch        *b_pos_y_det7;   //!
   TBranch        *b_pos_z_det7;   //!
   TBranch        *b_pos_x_target;   //!
   TBranch        *b_pos_y_target;   //!
   TBranch        *b_pos_z_target;   //!
   TBranch        *b_pos_x_dummy;   //!
   TBranch        *b_pos_y_dummy;   //!
   TBranch        *b_pos_z_dummy;   //!
   TBranch        *b_eloss_dummy;   //!
   TBranch        *b_inang_target;   //!
   TBranch        *b_Reaang_target;   //!
   TBranch        *b_ener_det0;   //!
   TBranch        *b_ener_det1;   //!
   TBranch        *b_ener_det2;   //!
   TBranch        *b_ener_det3;   //!
   TBranch        *b_ener_det4;   //!
   TBranch        *b_ener_det5;   //!
   TBranch        *b_ener_det6;   //!
   TBranch        *b_ener_det7;   //!
   TBranch        *b_rThetaCM0;   //!
   TBranch        *b_rThetaCM1;   //!
   TBranch        *b_rThetaCM2;   //!
   TBranch        *b_rThetaCM3;   //!
   TBranch        *b_rThetaCM4;   //!
   TBranch        *b_rThetaCM5;   //!
   TBranch        *b_rThetaCM6;   //!
   TBranch        *b_rThetaCM7;   //!
   TBranch        *b_t_sili0;   //!
   TBranch        *b_t_sili1;   //!
   TBranch        *b_t_sili2;   //!
   TBranch        *b_t_sili3;   //!
   TBranch        *b_t_sili4;   //!
   TBranch        *b_t_sili5;   //!
   TBranch        *b_t_sili6;   //!
   TBranch        *b_t_sili7;   //!
   TBranch        *b_px_dssd;   //!
   TBranch        *b_py_dssd;   //!
   TBranch        *b_pz_dssd;   //!
   TBranch        *b_Strip_Number0;   //!
   TBranch        *b_Strip_Number1;   //!
   TBranch        *b_Strip_Number2;   //!
   TBranch        *b_Strip_Number3;   //!
   TBranch        *b_Strip_Number4;   //!
   TBranch        *b_Strip_Number5;   //!
   TBranch        *b_Strip_Number6;   //!
   TBranch        *b_Strip_Number7;   //!
   TBranch        *b_truthPosx;   //!
   TBranch        *b_truthPosy;   //!
   TBranch        *b_truthPosz;   //!
   TBranch        *b_truthAngle_theta;   //!
   TBranch        *b_truthAngle_phi;   //!
   TBranch        *b_t_dssd;   //!
   TBranch        *b_t_dssd2;   //!
   TBranch        *b_Etot;   //!

   AnaClass3(TTree *tree=0);
   virtual ~AnaClass3();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnaClass3_cxx
AnaClass3::AnaClass3(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("SiTelescope",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("SiTelescope","");
      //chain->Add("tree_elas_12N_197Au.root/SiTelescope");
      chain->Add("tree_bu_12N_197Au.root/SiTelescope");
      //chain->Add("tree_elasXsec.root/SiTelescope");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

AnaClass3::~AnaClass3()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaClass3::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaClass3::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnaClass3::Init(TTree *tree)
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
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("E0", &E0, &b_E0);
   fChain->SetBranchAddress("E1", &E1, &b_E1);
   fChain->SetBranchAddress("E2", &E2, &b_E2);
   fChain->SetBranchAddress("E3", &E3, &b_E3);
   fChain->SetBranchAddress("E4", &E4, &b_E4);
   fChain->SetBranchAddress("E5", &E5, &b_E5);
   fChain->SetBranchAddress("E6", &E6, &b_E6);
   fChain->SetBranchAddress("E7", &E7, &b_E7);
   fChain->SetBranchAddress("EKin0", &EKin0, &b_EKin0);
   fChain->SetBranchAddress("EKin1", &EKin1, &b_EKin1);
   fChain->SetBranchAddress("EKin2", &EKin2, &b_EKin2);
   fChain->SetBranchAddress("EKin3", &EKin3, &b_EKin3);
   fChain->SetBranchAddress("EKin4", &EKin4, &b_EKin4);
   fChain->SetBranchAddress("EKin5", &EKin5, &b_EKin5);
   fChain->SetBranchAddress("EKin6", &EKin6, &b_EKin6);
   fChain->SetBranchAddress("EKin7", &EKin7, &b_EKin7);
   fChain->SetBranchAddress("pos_x_det0", &pos_x_det0, &b_pos_x_det0);
   fChain->SetBranchAddress("pos_y_det0", &pos_y_det0, &b_pos_y_det0);
   fChain->SetBranchAddress("pos_z_det0", &pos_z_det0, &b_pos_z_det0);
   fChain->SetBranchAddress("pos_x_det1", &pos_x_det1, &b_pos_x_det1);
   fChain->SetBranchAddress("pos_y_det1", &pos_y_det1, &b_pos_y_det1);
   fChain->SetBranchAddress("pos_z_det1", &pos_z_det1, &b_pos_z_det1);
   fChain->SetBranchAddress("pos_x_det2", &pos_x_det2, &b_pos_x_det2);
   fChain->SetBranchAddress("pos_y_det2", &pos_y_det2, &b_pos_y_det2);
   fChain->SetBranchAddress("pos_z_det2", &pos_z_det2, &b_pos_z_det2);
   fChain->SetBranchAddress("pos_x_det3", &pos_x_det3, &b_pos_x_det3);
   fChain->SetBranchAddress("pos_y_det3", &pos_y_det3, &b_pos_y_det3);
   fChain->SetBranchAddress("pos_z_det3", &pos_z_det3, &b_pos_z_det3);
   fChain->SetBranchAddress("pos_x_det4", &pos_x_det4, &b_pos_x_det4);
   fChain->SetBranchAddress("pos_y_det4", &pos_y_det4, &b_pos_y_det4);
   fChain->SetBranchAddress("pos_z_det4", &pos_z_det4, &b_pos_z_det4);
   fChain->SetBranchAddress("pos_x_det5", &pos_x_det5, &b_pos_x_det5);
   fChain->SetBranchAddress("pos_y_det5", &pos_y_det5, &b_pos_y_det5);
   fChain->SetBranchAddress("pos_z_det5", &pos_z_det5, &b_pos_z_det5);
   fChain->SetBranchAddress("pos_x_det6", &pos_x_det6, &b_pos_x_det6);
   fChain->SetBranchAddress("pos_y_det6", &pos_y_det6, &b_pos_y_det6);
   fChain->SetBranchAddress("pos_z_det6", &pos_z_det6, &b_pos_z_det6);
   fChain->SetBranchAddress("pos_x_det7", &pos_x_det7, &b_pos_x_det7);
   fChain->SetBranchAddress("pos_y_det7", &pos_y_det7, &b_pos_y_det7);
   fChain->SetBranchAddress("pos_z_det7", &pos_z_det7, &b_pos_z_det7);
   fChain->SetBranchAddress("pos_x_target", &pos_x_target, &b_pos_x_target);
   fChain->SetBranchAddress("pos_y_target", &pos_y_target, &b_pos_y_target);
   fChain->SetBranchAddress("pos_z_target", &pos_z_target, &b_pos_z_target);
   fChain->SetBranchAddress("pos_x_dummy", &pos_x_dummy, &b_pos_x_dummy);
   fChain->SetBranchAddress("pos_y_dummy", &pos_y_dummy, &b_pos_y_dummy);
   fChain->SetBranchAddress("pos_z_dummy", &pos_z_dummy, &b_pos_z_dummy);
   fChain->SetBranchAddress("eloss_dummy", &eloss_dummy, &b_eloss_dummy);
   fChain->SetBranchAddress("inang_target", &inang_target, &b_inang_target);
   fChain->SetBranchAddress("Reaang_target", &Reaang_target, &b_Reaang_target);
   fChain->SetBranchAddress("ener_det0", &ener_det0, &b_ener_det0);
   fChain->SetBranchAddress("ener_det1", &ener_det1, &b_ener_det1);
   fChain->SetBranchAddress("ener_det2", &ener_det2, &b_ener_det2);
   fChain->SetBranchAddress("ener_det3", &ener_det3, &b_ener_det3);
   fChain->SetBranchAddress("ener_det4", &ener_det4, &b_ener_det4);
   fChain->SetBranchAddress("ener_det5", &ener_det5, &b_ener_det5);
   fChain->SetBranchAddress("ener_det6", &ener_det6, &b_ener_det6);
   fChain->SetBranchAddress("ener_det7", &ener_det7, &b_ener_det7);
   fChain->SetBranchAddress("rThetaCM0", &rThetaCM0, &b_rThetaCM0);
   fChain->SetBranchAddress("rThetaCM1", &rThetaCM1, &b_rThetaCM1);
   fChain->SetBranchAddress("rThetaCM2", &rThetaCM2, &b_rThetaCM2);
   fChain->SetBranchAddress("rThetaCM3", &rThetaCM3, &b_rThetaCM3);
   fChain->SetBranchAddress("rThetaCM4", &rThetaCM4, &b_rThetaCM4);
   fChain->SetBranchAddress("rThetaCM5", &rThetaCM5, &b_rThetaCM5);
   fChain->SetBranchAddress("rThetaCM6", &rThetaCM6, &b_rThetaCM6);
   fChain->SetBranchAddress("rThetaCM7", &rThetaCM7, &b_rThetaCM7);
   fChain->SetBranchAddress("t_sili0", &t_sili0, &b_t_sili0);
   fChain->SetBranchAddress("t_sili1", &t_sili1, &b_t_sili1);
   fChain->SetBranchAddress("t_sili2", &t_sili2, &b_t_sili2);
   fChain->SetBranchAddress("t_sili3", &t_sili3, &b_t_sili3);
   fChain->SetBranchAddress("t_sili4", &t_sili4, &b_t_sili4);
   fChain->SetBranchAddress("t_sili5", &t_sili5, &b_t_sili5);
   fChain->SetBranchAddress("t_sili6", &t_sili6, &b_t_sili6);
   fChain->SetBranchAddress("t_sili7", &t_sili7, &b_t_sili7);
   fChain->SetBranchAddress("px_dssd", &px_dssd, &b_px_dssd);
   fChain->SetBranchAddress("py_dssd", &py_dssd, &b_py_dssd);
   fChain->SetBranchAddress("pz_dssd", &pz_dssd, &b_pz_dssd);
   fChain->SetBranchAddress("Strip_Number0", &Strip_Number0, &b_Strip_Number0);
   fChain->SetBranchAddress("Strip_Number1", &Strip_Number1, &b_Strip_Number1);
   fChain->SetBranchAddress("Strip_Number2", &Strip_Number2, &b_Strip_Number2);
   fChain->SetBranchAddress("Strip_Number3", &Strip_Number3, &b_Strip_Number3);
   fChain->SetBranchAddress("Strip_Number4", &Strip_Number4, &b_Strip_Number4);
   fChain->SetBranchAddress("Strip_Number5", &Strip_Number5, &b_Strip_Number5);
   fChain->SetBranchAddress("Strip_Number6", &Strip_Number6, &b_Strip_Number6);
   fChain->SetBranchAddress("Strip_Number7", &Strip_Number7, &b_Strip_Number7);
   fChain->SetBranchAddress("truthPosx", &truthPosx, &b_truthPosx);
   fChain->SetBranchAddress("truthPosy", &truthPosy, &b_truthPosy);
   fChain->SetBranchAddress("truthPosz", &truthPosz, &b_truthPosz);
   fChain->SetBranchAddress("truthAngle_theta", &truthAngle_theta, &b_truthAngle_theta);
   fChain->SetBranchAddress("truthAngle_phi", &truthAngle_phi, &b_truthAngle_phi);
   fChain->SetBranchAddress("t_dssd", &t_dssd, &b_t_dssd);
   fChain->SetBranchAddress("t_dssd2", &t_dssd2, &b_t_dssd2);
   fChain->SetBranchAddress("Etot", &Etot, &b_Etot);
   Notify();
}

Bool_t AnaClass3::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaClass3::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaClass3::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnaClass3_cxx
