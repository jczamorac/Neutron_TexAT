{
TChain chain("SiTelescope") ;
chain.Add("tree.root");
chain.GetListOfFiles()->Print();
chain.MakeClass("AnaClass3");
}
