#include "interface/HiggsCrossSectionReader.hh"
#include "interface/MiniTreeFitter1D.hh"

#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TRint.h>

#include <iomanip>
#include <boost/program_options.hpp>

using namespace std;
using namespace RooFit;


int main( int nargc, char **argv ) {

  namespace po = boost::program_options;
  string config_file;

  int IJOB, iCosThetaStarBin;
  string dir;
  po::options_description config("Configuration");
  config.add_options()
    ("help,h"   ,"help")
    ("iCSbin,b" ,po::value<int>(&iCosThetaStarBin)->default_value(0),"cosThetaStar bin" )
    ("ijob,j"   ,po::value<int>(&IJOB)            ->default_value(0),"job number")
    ("inDir,d"  ,po::value<string>(&dir)          ->default_value("/Users/fcouderc/CMSWork/h2gg/diphoton2012_mvaSel_cms53x_v10/"),"directory with input")
    ;

  po::variables_map vm;
  po::store(po::command_line_parser(nargc, argv).
	    options(config).run(), vm);
  po::notify(vm);

  if( vm.count("help") ) {
    cout << config << endl;
    cout << "Usage: ./bin/OptiMvaPhoID [options]" << endl;
    return 1;
 }

  TRint interactive("rootInteractive",(int*)0,0);
  
  vector<TCut> cut_cos;
  cut_cos.push_back("abs(cThetaStar_CS)>=0.0 && abs(cThetaStar_CS)<0.2");
  cut_cos.push_back("abs(cThetaStar_CS)>=0.2 && abs(cThetaStar_CS)<0.4");
  cut_cos.push_back("abs(cThetaStar_CS)>=0.4 && abs(cThetaStar_CS)<0.6");
  cut_cos.push_back("abs(cThetaStar_CS)>=0.6 && abs(cThetaStar_CS)<0.8");
  cut_cos.push_back("abs(cThetaStar_CS)>=0.8 && abs(cThetaStar_CS)<1.0");

  int it10 = 4; // <=> -0.8

  cout << " ======== I will do jobs " << endl;
  int ijob = -1;
  for( int it1 = it10 ; it1 < 39; it1++ ) 
  for( int it2 = it1+1; it2 < 40; it2++ )  {
    ijob++;
    if( ijob%50 != IJOB ) continue;
    cout << "   job = " << ijob << endl;
  }
  cout << " # of job = " << ijob << endl;

  ijob = -1;  
  for( int it1 = it10 ; it1 < 39; it1++ ) 
  for( int it2 = it1+1; it2 < 40; it2++ )  {
    ijob++;
    if( ijob%50 != IJOB ) continue;


    cout << "=============== it1: " << it1 << " ; it2: " << it2 << " ================" << endl; 

    float cut1 = -0.2 + it1*0.05;
    float cut2 = -0.2 + it2*0.05;
    char cutmva1[1000], cutmva2[1000];
    char directory[1000];
    sprintf(cutmva1,"diphoMva>=%1.3f&&diphoMva<%1.3f",cut1,cut2);
    sprintf(cutmva2,"diphoMva>=%1.3f&&diphoMva<99999",cut2);
    sprintf(directory,"results/optimizeSpin/cThetaStar:%d/cut1:%1.3f_cut2:%1.3f/",iCosThetaStarBin,cut1,cut2);
    
    vector<TCut> cut_category;
    cut_category.push_back( TCut(cutmva1) && cut_cos[iCosThetaStarBin] ); 
    cut_category.push_back( TCut(cutmva2) && cut_cos[iCosThetaStarBin] );
    
    //----  configure fitter (define output dir, categories, add variables...)
    string hlFactoryCard = "etc/workspaceConfig/mva_sm_model_cond.rs";
    MiniTreeFitter1D fitter(hlFactoryCard);
    fitter.setPlotDirectory(directory);
    fitter.setMainCut( TCut("mass>100&&mass<180") );
    fitter.setCategories( cut_category );
    RooRealVar *diphoMva  = new RooRealVar( "diphoMva"    ,"",-0.8,1.0); fitter.addVariable( diphoMva );
    RooRealVar *cosThetaS = new RooRealVar("cThetaStar_CS","",-1.5,1.5); fitter.addVariable(cosThetaS);
    
    bool addSignal = 1;
    bool addData   = 1;
    bool dofits    = 1;
    
    
    if( addSignal ) {
      float xsec = 20.00*1000; //fb
      float br = 0.00228;
      float lumi = 12.2; //fb-1
      
      vector<string> sFiles1, sNames1;
      vector<float>  sXsec1;
      sFiles1.push_back( dir + "/mc/job_hgg_gg0odd.root_0.root" );  sNames1.push_back( "spin 0+" ); sXsec1.push_back(xsec);
      fitter.addSigSamples(sFiles1,sXsec1,sNames1,br,lumi);
      
	if( dofits ) {
	  fitter.modelSignal(125,"Spin0_cThetaStar" + itostr(iCosThetaStarBin),0);
	  fitter.makeSignalWorkspace();
	}
    }
    
    if( addData ) {
      string dataset = dir + "/data/MiniTreeData.root";
      fitter.addData(dataset);
      if( dofits ) { 	
	fitter.modelBackground(125.,"_cThetaStar"+itostr(iCosThetaStarBin) );
	fitter.makeBackgroundWorkspace();
      }
    }
    
    if( addSignal && dofits )
      fitter.simultaneousFitOnlyOneSig(0);
    
    if( addData && addSignal && dofits ) {
      string datacard = "spin_mh125_cThetaStarBin" + itostr(iCosThetaStarBin) + ".datacard.txt";
      fitter.createDataCard(datacard);    
    }
  }

  return 1;
}



