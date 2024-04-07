#define AnaClass3_cxx
#include "AnaClass3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnaClass3::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L AnaClass2.C
  //      root> AnaClass2 t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
     if (fChain == 0) return;

     Long64_t nentries = fChain->GetEntriesFast();


     //TFile *file = new TFile("tamu_sim_elas_Xsec.root", "recreate");
     TFile *file = new TFile("tamu_sim_bu_Xsec.root", "recreate");
     file->mkdir("Det1");
     file->mkdir("Det2");
     file->mkdir("Det3");

     int hbin = 3*64;

     TH1F *EKin_45 = new TH1F("Ekin_45", "", 360, 0, 90);
     TH2F *Kplot_45 = new TH2F("Kplot_45", "", hbin, 0, 180, 200,2,90);
     TH1F *EKin_90 = new TH1F("Ekin_90", "", 360, 0, 90);
     TH2F *Kplot_90 = new TH2F("Kplot_90", "", hbin, 0, 180, 200,2,90);
     TH1F *EKin_130 = new TH1F("Ekin_130", "", 90, 0, 90);
     TH2F *Kplot_130 = new TH2F("Kplot_130", "", hbin, 0, 180, 200,2,90);
     /*
     TH1F *Thcmout = new TH1F("Thcmout", "", 200, 120, 180);
     TH1F *Kinout = new TH1F("Kinout", "", 500, 0, 20);

     TH2F *Bipara0 = new TH2F("Bipara0", "", 60, 0, 60, 500,0,20);
     TH2F *Bipara1 = new TH2F("Bipara1", "", 60, 0, 60, 500,0,20);
     TH2F *Bipara2 = new TH2F("Bipara2", "", 60, 0, 60, 500,0,20);
     TH2F *Bipara3 = new TH2F("Bipara3", "", 60, 0, 60, 500,0,20);
     TH2F *BiparaTB = new TH2F("BiparaTB", "", 60, 0, 60, 500,0,20);
     TH2F *Bipara4 = new TH2F("Bipara4", "", 60, 0, 60, 500,0,20);
     TH2F *Bipara5 = new TH2F("Bipara5", "", 60, 0, 60, 500,0,20);
     TH2F *Bipara6 = new TH2F("Bipara6", "", 60, 0, 60, 500,0,20);
     TH2F *Bipara7 = new TH2F("Bipara7", "", 60, 0, 60, 500,0,20);
     TH2F *BiparaTF = new TH2F("BiparaTF", "", 60, 0, 60, 500,0,20);

     TH2F *ThCMStr = new TH2F("ThCMStr", "", 60, 0, 60, 180,0,0);
     */

     TH1F *ESpec_det1[9];
     TH1F *ESpec_det2[9];
     TH1F *ESpec_det3[9];

     for(int i=0; i<9;i++){
	      float ang1 = i*5+22.5;
        float ang2 = i*5+67.5;
        float ang3 = i*5+97.5;
	      ESpec_det1[i] =  new TH1F(Form("det1_th_%.1f",ang1),"",90,0,90);
	      ESpec_det2[i] =  new TH1F(Form("det2_th_%.1f",ang2),"",90,0,90);
        ESpec_det3[i] =  new TH1F(Form("det3_th_%.1f",ang3),"",90,0,90);
      }



     Long64_t nbytes = 0, nb = 0;
     for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        int nbinx = 64;
        int nbiny = 64;
        double step = 100./nbinx;
        double xmin = -50;
        double ymin = -50;


        double ang_det1 = -45;
        double ang_det2 = -90;
        double ang_det3 = -120;
        double rot_det = -45;


        double z2 = 153;
        double z1 = 153;
        double z0 = 153;
        double x,y, z, r, theta, phi;
        double xrot, yrot;
        //ouble factor = 1;
        double factor = 1/128.0;



        if (E0>0){

          int XBin = TMath::Floor( (pos_x_det0 - xmin) / step );
          int YBin = TMath::Floor( (pos_y_det0 - ymin) / step );
          double Xnew = xmin + step*(XBin+0.5);
          double Ynew = ymin + step*(YBin+0.5);
          r = sqrt(z0*z0 + Xnew*Xnew + Ynew*Ynew);

          //rotation diamont config
          xrot = Xnew*cos((rot_det)*3.141562/180.) - (Ynew)*sin((rot_det)*3.141562/180.);
          yrot = Xnew*sin((rot_det)*3.141562/180.) + (Ynew)*cos((rot_det)*3.141562/180.);


          x = xrot*cos((ang_det1)*3.141562/180.) + z0*sin((ang_det1)*3.141562/180.);
          y = yrot;
          z = -(xrot)*sin((ang_det1)*3.141562/180.) + z0*cos((ang_det1)*3.141562/180.);

          theta = acos(z/r)*180/3.1415;
          Kplot_45->Fill(theta,E0, factor);

          if(theta>=40 && theta<=45)EKin_45->Fill(E0, factor);
          //cout<<pos_x_det0<<"  "<<Xnew<<"  "<<pos_y_det0<<"  "<<Ynew<<"  "<<theta<<endl;
          if(theta>20 && theta<65){
            int de1Bin = TMath::Floor( (theta - 20.0) / 5 );
            ESpec_det1[de1Bin]->Fill(E0, factor);
          }
        }

        if (E1>0){
          //EKin_90->Fill(E1);
          int XBin = TMath::Floor( (pos_x_det1 - xmin) / step );
          int YBin = TMath::Floor( (pos_y_det1 - ymin) / step );
          double Xnew = xmin + step*(XBin+0.5);
          double Ynew = ymin + step*(YBin+0.5);
          r = sqrt(z1*z1 + Xnew*Xnew + Ynew*Ynew);

          //rotation diamont config
          xrot = Xnew*cos((rot_det)*3.141562/180.) - (Ynew)*sin((rot_det)*3.141562/180.);
          yrot = Xnew*sin((rot_det)*3.141562/180.) + (Ynew)*cos((rot_det)*3.141562/180.);


          x = xrot*cos((ang_det2)*3.141562/180.) + z1*sin((ang_det2)*3.141562/180.);
          y = yrot;
          z = -(xrot)*sin((ang_det2)*3.141562/180.) + z1*cos((ang_det2)*3.141562/180.);

          theta = acos(z/r)*180/3.1415;
          Kplot_90->Fill(theta,E1, factor);
          if(theta>=80 && theta<=85)EKin_90->Fill(E1, factor);
          //cout<<x<<"  "<<y<<"  "<<z<<"  "<<sqrt(x*x+y*y+z*z)<<"  "<<theta<<endl;
          if(theta>65 && theta<105){
            int de2Bin = TMath::Floor( (theta - 65.0) / 5 );
            ESpec_det2[de2Bin]->Fill(E1, factor);
          }
        }

        if (E2>0){
          //EKin_130->Fill(E2);
          int XBin = TMath::Floor( (pos_x_det2 - xmin) / step );
          int YBin = TMath::Floor( (pos_y_det2 - ymin) / step );
          double Xnew = xmin + step*(XBin+0.5);
          double Ynew = ymin + step*(YBin+0.5);
          r = sqrt(z2*z2 + Xnew*Xnew + Ynew*Ynew);


          //rotation diamont config
          xrot = Xnew*cos((rot_det)*3.141562/180.) - (Ynew)*sin((rot_det)*3.141562/180.);
          yrot = Xnew*sin((rot_det)*3.141562/180.) + (Ynew)*cos((rot_det)*3.141562/180.);


          x = xrot*cos((ang_det3)*3.141562/180.) + z2*sin((ang_det3)*3.141562/180.);
          y = yrot;
          z = -(xrot)*sin((ang_det3)*3.141562/180.) + z2*cos((ang_det3)*3.141562/180.);

          theta = acos(z/r)*180/3.1415;
          Kplot_130->Fill(theta,E2, factor);
          if(theta>=110 && theta<=115)EKin_130->Fill(E2, factor);
          //cout<<x<<"  "<<y<<"  "<<z<<"  "<<sqrt(x*x+y*y+z*z)<<"  "<<theta<<endl;
          if(theta>95 && theta<125){
            int de3Bin = TMath::Floor( (theta - 95.0) / 5 );
            ESpec_det3[de3Bin]->Fill(E2, factor);
          }
        }


      }

     EKin_45->Write();
     Kplot_45->Write();
     EKin_90->Write();
     Kplot_90->Write();
     EKin_130->Write();
     Kplot_130->Write();
     for(int i=0; i<9;i++){
	      file->cd("Det1");
	      ESpec_det1[i]->Write();
        file->cd("Det2");
	      ESpec_det2[i]->Write();
        file->cd("Det3");
	      ESpec_det3[i]->Write();
      }



     file->Close();
}
