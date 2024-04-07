{
gStyle->SetTitleBorderSize(0);
  gStyle->SetFrameFillColor(kWhite); // needed when using root v5
  gStyle->SetPadColor(kWhite);
  gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping
  //gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9
  gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plott
  // For publishing:
  gStyle->SetLineWidth(2.);
  gStyle->SetTextSize(1.1);
  gStyle->SetLabelSize(0.06,"xy");
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.12);


gROOT->ProcessLine(".L $ROOTSYS/lib/libPhysics.so");
gROOT->ProcessLine(".L AnaClass3.C+");
gROOT->ProcessLine("AnaClass3 t");
}
