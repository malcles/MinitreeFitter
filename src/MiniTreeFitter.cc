#include "HiggsCrossSectionReader.hh"
#include "MiniTreeFitter1D.hh"


#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TRint.h>

#include <iomanip>
#include <boost/program_options.hpp>
using namespace std;



int main( int nargc, char **argv ) {
    
  namespace po = boost::program_options;
  string config_file;
  string dirAFS;
  float mMin, mMax, mh;
  int  bkgModel;
  int  cat;
  bool addSig, addBkg, simFit;
  bool doFits = true;
  bool rootInteractive(true);
  po::options_description config("Configuration");
  config.add_options()
    ("help,h"   ,"help")
    ("mh"  , po::value<float>(&mh)       ->default_value(125),"Higgs Mass" )
    ("bkg" , po::value<int>(&bkgModel)   ->default_value(0)  ,"bkg model")
    ("mMin",po::value<float>(&mMin)      ->default_value(100.),"minimum mass cut")
    ("mMax",po::value<float>(&mMax)      ->default_value(180.),"maximum mass cut")
    //("inDir,d",po::value<string>(&dirAFS)->default_value("../diphoton2012_mvaSel_cms53x_v10/"),"directory with input")
    // JM my data dir:
    ("inDir,d",po::value<string>(&SWdirAFS)->default_value("/afs/cern.ch/user/m/malcles/MonScratch/private/MiniTreeFitter/data/CiC/"),"directory with input") 
    ("addSig", po::value<bool>(&addSig)  ->default_value(true ),"add signal" )
    ("addBkg", po::value<bool>(&addBkg)  ->default_value(true ),"add background" )
    ("simFit", po::value<bool>(&simFit)  ->default_value(false),"simultaneous fit" )
    ("cat"   , po::value<int>(&cat)      ->default_value(-1),"category" )
    ("doFits", po::value<bool>(&doFits)      ->default_value(1),"do fits" )
    ("interactive,i", po::value<bool>(&rootInteractive)->default_value(1),"root interactive" )

    ;
  
  po::variables_map vm;SW
  po::store(po::command_line_parser(nargc, argv).
	    options(config).run(), vm);
  po::notify(vm);
  
  if( vm.count("help") ) {
    cout << config << endl;
    cout << "Usage: ./bin/MiniTreeFitter [options]" << endl;
    return 1;
  }
  
  TRint *interactive = 0;
  if( rootInteractive ) interactive = new TRint("rootInteractive",(int*)0,0);

  string categorisation = "notDefined";
  if( dirAFS.find( "mva" ) != string::npos ||  dirAFS.find( "MVA" ) != string::npos ) categorisation = "mva";
  if( dirAFS.find( "cic" ) != string::npos ||  dirAFS.find( "CiC" ) != string::npos ||
      dirAFS.find( "CIC" ) != string::npos )  categorisation = "cicpf";
  
  if( categorisation == "notDefined" ) {
    cout << "  categorisation not defined: directory name should contain mva or cic or cicpf" << endl;
    return 1;
  }



  // ----  Define the categories ---- //
  vector<TCut> smCategories;
  vector<TString> smCatNames;
  vector<int>  polOrder;
  if( categorisation == "mva" ) {
    smCategories.push_back( "catMva == 0" ); polOrder.push_back( 5 );
    smCategories.push_back( "catMva == 1" ); polOrder.push_back( 5 );
    smCategories.push_back( "catMva == 2" ); polOrder.push_back( 5 );
    smCategories.push_back( "catMva == 3" ); polOrder.push_back( 5 );   
  } else if( categorisation == "cicpf" ) {
    //smCategories.push_back( "tagCat == 8 && abs(dEtaJJ) > 3" ); polOrder.push_back( 3 );
    //smCategories.push_back( "tagCat == 8 && abs(dEtaJJ) < 3" ); polOrder.push_back( 5 );   

    // JM: try more categories with my ntuples:
    //==========================================
    smCategories.push_back( "tagCat < 0 && untagCat==0 " ); polOrder.push_back( 4 ); smCatNames.push_back("untag_0");  // untag 0
    smCategories.push_back( "tagCat < 0 && untagCat==1 " ); polOrder.push_back( 5 ); smCatNames.push_back("untag_1");  // untag 1
    smCategories.push_back( "tagCat < 0 && untagCat==2 " ); polOrder.push_back( 4 ); smCatNames.push_back("untag_2");  // untag 2
    smCategories.push_back( "tagCat < 0 && untagCat==3 " ); polOrder.push_back( 5 ); smCatNames.push_back("untag_3");  // untag 3
    smCategories.push_back( "tagCat < 0 && untagCat==4 " ); polOrder.push_back( 4 ); smCatNames.push_back("untag_4");  // untag 4
    smCategories.push_back( "tagCat < 0 && untagCat==5 " ); polOrder.push_back( 5 ); smCatNames.push_back("untag_5"); // untag 5
    smCategories.push_back( "tagCat < 0 && untagCat==6 " ); polOrder.push_back( 4 ); smCatNames.push_back("untag_6"); // untag 6
    smCategories.push_back( "tagCat < 0 && untagCat==7 " ); polOrder.push_back( 5 ); smCatNames.push_back("untag_7"); // untag 7
    smCategories.push_back( "tagCat == 8 " ); polOrder.push_back( 3 ); smCatNames.push_back("vbf_tig"); // vbf 
    smCategories.push_back( "tagCat == 9 " ); polOrder.push_back( 3 ); smCatNames.push_back("vbf_loo"); // vbf
    smCategories.push_back( "tagCat == 10 " ); polOrder.push_back( 2 ); smCatNames.push_back("vhl_tig"); // vh
    smCategories.push_back( "tagCat == 11 " ); polOrder.push_back( 4 ); smCatNames.push_back("vhl_loo"); // vh
    smCategories.push_back( "tagCat == 12 " ); polOrder.push_back( 3 ); smCatNames.push_back("vh_met "); // vh
    smCategories.push_back( "tagCat == 13 " ); polOrder.push_back( 2 ); smCatNames.push_back("tth_lep"); // tth
    smCategories.push_back( "tagCat == 14" ); polOrder.push_back( 2 ); smCatNames.push_back("tth_had"); // tth
    smCategories.push_back( "tagCat == 15 " ); polOrder.push_back( 3 ); smCatNames.push_back("vh_had "); // vhhad
  }

  vector<int> polOrderCuts;
  vector<TCut> smCatCuts;
  if( cat >=  0 && cat < int(smCategories.size()) ) {
    smCatCuts.push_back( smCategories[cat] );
    polOrderCuts.push_back(polOrder[cat]);
  }
  else  {
    smCatCuts = smCategories;
    polOrderCuts = polOrder;
  }

  // ----  prepare the fitter ---- //
  string hlFactoryCard = "etc/workspaceConfig/mva_sm_model_cond.rs";
  MiniTreeFitter1D fitter(hlFactoryCard);
  fitter.setPlotDirectory("workspace_SM/" + categorisation + "/");
  // fitter.setMainCut( basicCut );

  string mName = "mass";
  //  RooRealVar *catMva   = new RooRealVar( "catMva"  ,"",-999,10.0); fitter.addVariable( catMva   );
  //  RooRealVar *diphoMva = new RooRealVar( "diphoMva","",-3  ,10.0); fitter.addVariable( diphoMva );
  RooRealVar *catBase = new RooRealVar( "tagCat" ,"",-999,30.0); fitter.addVariable( catBase  );
  RooRealVar *catBaseUntag = new RooRealVar( "untagCat" ,"",-999,30.0); fitter.addVariable( catBaseUntag  );
  RooRealVar *dEtaJJ  = new RooRealVar( "dEtaJJ" ,"",-999,10.0); fitter.addVariable( dEtaJJ  );
  RooRealVar *mass    = new RooRealVar( mName.c_str(),"",mMin,mMax); fitter.addVariable( mass     );
  fitter.setMassVarName( mName );
  fitter.setMassVarSet(true);

  fitter.setCategories( smCatCuts );
  fitter.setCategoriesNames( smCatNames );
  fitter.setMassMin(mMin);
  fitter.setMassMax(mMax);

  string escaleSyst = "";  
  string categorySuffix = "_cat" + itostr(cat);
  string dirData = dirAFS + "/data/";
  string dirMC   = dirAFS + "/mc/" + escaleSyst;

  if( addSig ) {
    SMHiggsCrossSection HiggsXS; HiggsXS.is8TeV();
    float xsec_ggh = 1000*HiggsXS.HiggsSMxsec_ggh(mh); //fb
    float xsec_vbf = 1000*HiggsXS.HiggsSMxsec_vbf(mh); //fb
    float xsec_vh  = 1000*(HiggsXS.HiggsSMxsec_wh(mh) + HiggsXS.HiggsSMxsec_zh(mh)); //fb
    float xsec_tth = 1000*HiggsXS.HiggsSMxsec_tth(mh); //fb
    float br = HiggsXS.HiggsBR(mh);
    float lumi = 19.5; //fb-1
    cout << " br = " << br << " - xsec_ggh = " << xsec_ggh << endl;
    vector<string> sFiles1, sNames1;
    vector<float>  sXsec1;
    //    sFiles1.push_back( dirMC + "/mc/job_hgg_gg0odd.root_0.root" );  sNames1.push_back( "spin 0+" ); sXsec1.push_back(20*1000);
    sFiles1.push_back( dirMC + TString::Format("job_summer12_ggH_%3.0f.root_0.root",mh).Data()  );  sNames1.push_back( "ggH"); sXsec1.push_back(xsec_ggh);
    sFiles1.push_back( dirMC + TString::Format("job_summer12_VBF_%3.0f.root_0.root",mh).Data()  );  sNames1.push_back( "VBF"); sXsec1.push_back(xsec_vbf);
    sFiles1.push_back( dirMC + TString::Format("job_summer12_WH_ZH_%3.0f.root_0.root",mh).Data());  sNames1.push_back( "VH" ); sXsec1.push_back(xsec_vh );
    sFiles1.push_back( dirMC + TString::Format("job_summer12_TTH_%3.0f.root_0.root",mh).Data()  );  sNames1.push_back( "ttH"); sXsec1.push_back(xsec_tth); 

    
    fitter.addSigSamples(sFiles1,sXsec1,sNames1,br,lumi);
    if( doFits) {    
      fitter.modelSignal(125,categorySuffix);
      fitter.makeSignalWorkspace();
    }
  }
  
  float signi = -1;
  vector<double> signis;
  if( addBkg ) {

    //    string dataset = dirData + "MiniTreeData_2012abcd.root";
    string dataset = dirData + "data_8TeV_skimMVA_runABCD.root";

    fitter.unblind();
    fitter.addData(dataset);
    if( doFits ) {
      if( bkgModel == 0 ){
	signis = fitter.modelBackground(125.,categorySuffix); 
	fitter.makeBackgroundWorkspace();
      }
      else if (bkgModel == 1) { 
	fitter.setPolynomialOrder(  polOrderCuts );
	fitter.modelBackground(   125.,categorySuffix );
	fitter.modelBackgroundPol(125.,categorySuffix );
	fitter.modelBackgroundExp(125.,categorySuffix );
	// fitter.makeBackgroundWorkspace("Pol");
	  // fitter.makeBackgroundWorkspace();
	fitter.dumpBkgFitParam();
	if( simFit ) {
	  fitter.backgroundFitResut(cat,0);
	}
      }
    }
  }
  //  return;
  if( addSig && addBkg && simFit ) fitter.simultaneousFitOnlyOneSig(0);

  if( addSig && addBkg && doFits ) {
    string datacard = "MVAFact_SM_mh125_cat" + itostr(cat) + ".datacard.txt";
    fitter.createDataCard(datacard);    
  }
  signi = 0;
  for( unsigned ic = 0 ;ic < signis.size(); ic++ ) {
    cout << " signi cat[" << ic << "] = " << signis[ic] << endl;
    signi += signis[ic]*signis[ic];
  }
  cout << " SigniTot: " << sqrt(signi) << endl;

  if( interactive ) {
    interactive->Run();
    delete interactive;
  }
  
  return 0;
}

