#include <TString.h>
#include <TTree.h>
#include <TFile.h>

#include <iostream>
#include <string>
using namespace std;


void ChangeWeight(string rootfile) {

  float scale1 = 10;
  float scale2 = 0.1;
  
  if( rootfile.find("Spin2") != string::npos ) {
    scale1 = 0.1;
    scale2 = 10;
  }

  vector<string> brToActivate;
  brToActivate.push_back("mass");
  brToActivate.push_back("maxSCEta");
  brToActivate.push_back("minR9");
  brToActivate.push_back("catBase");
  brToActivate.push_back("catMva");
  brToActivate.push_back("diphoMva");
  brToActivate.push_back("massResoEng");
  brToActivate.push_back("massDefVtx");
  brToActivate.push_back("cThetaStar_CS");

  /// input tree
  TTree *tin = (TTree*) TFile::Open(rootfile.c_str(),"read")->Get("HToGG");
  tin->SetBranchStatus("*",0);
  for( unsigned ibr = 0 ; ibr < brToActivate.size(); ibr++ )   tin->SetBranchStatus(brToActivate[ibr].c_str(),1);

  /// output tree
  TFile *f = new TFile( string(rootfile + ".newModel.root").c_str(),"recreate");
  TTree *tout = tin->CloneTree();

  /// set branches 

  float wei, cTS;
  tin ->SetBranchStatus( "wei" , 1);
  tin ->SetBranchAddress("wei" , &wei);
  tin ->SetBranchAddress("cThetaStar_CS" ,&cTS);
  TBranch *brWei = tout->Branch("wei" ,&wei ,"wei/F" ); 

  int nentries = tin->GetEntries();
  cout << " -- Doing Tree: " << rootfile << " #evts = " << nentries << endl;
  for( int ientry = 0; ientry < nentries; ientry++ ) {
    tin->GetEntry(ientry);
    if( fabs( cTS ) < 0.45 ) wei *= scale1;
    else                     wei *= scale2;
    
    brWei->Fill();
  }
  
  tout->Write("HToGG",TObject::kOverwrite);
  f->Close();
}
