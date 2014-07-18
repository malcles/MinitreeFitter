#include "interface/SpinFitterConfig.hh"
#include "interface/SetupReader.hh"

#include <TCut.h>
#include <TString.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

vector<TCut> SpinFitterConfig::getCosThetaCuts(void) const {
  vector<TCut> cuts;

  for( unsigned ib = 0; ib < _cosThetaStarBoundaries.size()-1; ++ib ) {
    cuts.push_back( TCut(TString::Format("abs(cThetaStar_CS)>=%f&&abs(cThetaStar_CS)<%f",
					 _cosThetaStarBoundaries[ib],
					 _cosThetaStarBoundaries[ib+1])) );
  }
  
  return cuts;
}

vector<vector<TCut> > SpinFitterConfig::getDiphoDiscriCuts(void) const {
  vector<TCut> cosThetaCuts = getCosThetaCuts();
  vector<vector<TCut> > cuts;
  
  for( unsigned ib = 0; ib < _cosThetaStarBoundaries.size()-1; ++ib ) {
    vector<TCut> cuttmp;
    for( unsigned id = 0; id < _cosDiphoDiscriBoundaries[ib].size()-1; ++id ) {
      TCut c = TCut(TString::Format("diphoMva>=%f&&diphoMva<%f",
				    _cosDiphoDiscriBoundaries[ib][id],
				    _cosDiphoDiscriBoundaries[ib][id+1]));
      cuttmp.push_back( c && cosThetaCuts[ib] );
    }
    cuts.push_back( cuttmp );
  }
  return cuts;
}

SpinFitterConfig::SpinFitterConfig(string configname) { _config = configname; setup(); }
void SpinFitterConfig::setup() {
  SetupReader reader(_config);
  char varName[1000];

  bool found(true);
  unsigned ib = 0;
  while(found&&ib<1000) {
    sprintf(varName,"cosThetaBound%d",ib);
    float boundary = reader.getFloat( varName, found );
    if( found ) _cosThetaStarBoundaries.push_back( boundary );
    ib++;
  }
  
  unsigned id = 0;
  for( ib = 0; ib < _cosThetaStarBoundaries.size()-1; ++ib ) {
    vector<float> btmp; 
    found = true; id = 0;
    while(found&&id<1000) {
      sprintf(varName,"cT%d_mvaBound%d",ib,id);
      float boundary = reader.getFloat( varName, found );
      if( found ) btmp.push_back( boundary );
      id++;
    }
    _cosDiphoDiscriBoundaries.push_back( btmp );
  }
  
  vector<vector<TCut> > cuts = getDiphoDiscriCuts();
  cout << "========================= SpinFitter Config =======================" << endl;
  for(  ib = 0; ib < cuts.size(); ++ib ) {
    cout << "** cosThetaStar bin: "<< ib << "--------------" << endl;
    for(  id = 0; id < cuts[ib].size(); ++id ) {
      cout << cuts[ib][id] << endl;
    }
  }
  cout << "=========================+++++++++++++++++++=======================" << endl;

}

inline float SpinFitterConfig::AccTot(string rootfile, string treename ) const {
  TTree *t =0; t = (TTree*) TFile::Open(rootfile.c_str())->Get(treename.c_str());
  float ntot = 0;
  if( t ) {
    vector<vector<TCut> > cuts = getDiphoDiscriCuts();
    for( unsigned ic1 = 0 ; ic1 < cuts.size()     ; ++ic1 )
    for( unsigned ic2 = 0 ; ic2 < cuts[ic1].size(); ++ic2 ) {
      TH1F h("htmpMass","htmpMass",80,100,180);
      t->Draw("mass>>htmpMass","wei"*cuts[ic1][ic2],"goff");
      ntot += h.Integral();
      cout << "AccCat: " << h.Integral()/(1e5) << "  <-->  " << cuts[ic1][ic2] << endl;
    }
  }
  return ntot/float(1e5);
}

