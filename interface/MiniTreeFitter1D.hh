#ifndef MiniTreeFitter1D__hh__
#define MiniTreeFitter1D__hh__

#include <TCut.h>
#include <TSystem.h>

#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/HLFactory.h"
#include "RooStats/RooStatsUtils.h"
#include "RooSimultaneous.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

inline std::string itostr(int i ) {
  std::stringstream out;
  out << i;
  return out.str();
}

inline std::string ftostr(float f ) {
  std::stringstream out;
  out << f;
  return out.str();
}


class MiniTreeFitter1D {
public:
  inline MiniTreeFitter1D(std::string hlFactory);
  inline ~MiniTreeFitter1D(void);
  inline void setMainCut( TCut cut ) { _mainCut = cut;}
  inline void setCategories( std::vector<TCut> cuts ) { _categories = cuts; }
  inline void setCategoriesNames( std::vector<TString> names ) { _categoriesNames = names; } // JM
  inline void setPolynomialOrder(std::vector<int> polOrder) {_polynomialOrder = polOrder; }
  inline void setPlotDirectory( std::string dir );
  inline void addVariable( RooRealVar *v ) {  _variables.add(*v); }
  inline void addSyst( bool addsyst = true) { _addSyst = addsyst; } 

  //-- simultaneousFit
  void simultaneousFit(int sHyp);
  void simultaneousFitOnlyOneSig(int sHyp = 0 );

  //-- unset the blinding flag
  inline void unblind(void) { _unblind = true; }
  
  //--- add acceptance correction for alternate signal hypothesis
  inline void setAccCorrForAltSig( float corr ) {  _altSigAccCorr = corr; }

  //-- add signal samples: can add different signal hypothesis
  void addSigSamples(std::vector<std::string> samples    , std::vector<float> xsec, 
		     std::vector<std::string> signalNames, float br = 1, float lumi = 1);
  //-- fit the different signal hypothesis (by default last hypothesis added is fitted
  //-- 2 first params are just indicative (for saving xchecks plots and naming workspaces)
  void modelSignal( float testMass = 0, std::string sigHypo = "", int sHyp = -1 ) ;

  ///-- add data, this will be used for bckg shape fit
  ///-- final fit is performed in a reduced range around test mass [-20,20]
  void addData(std::string dataset);
  std::vector<double> modelBackground( float testMass = 0, std::string suffixName = "" ) ;
  double modelBackgroundPol( float testMass = 0, std::string suffixName = "" ) ;
  double modelBackgroundExp( float testMass = 0, std::string suffixName = "" ) ;
  void   makeBackgroundWorkspace( std::string bkg = "Pow" );
  void backgroundFitResut(int icat = -1, unsigned iBkgModelRef = 0 );

  void   dumpBkgFitParam(void);

  ///-- create datacard / workspaces
  void createDataCard( std::string cardname, bool useAltSig = false );
  void makeSignalWorkspace(void);

  inline void setMassMin(float m) { _fitMassMin = m; };
  inline void setMassMax(float m) { _fitMassMax = m; };
  inline std::string massVarName(void)    { return _massVarName; }
  inline void  setMassVarName( std::string mName ) { _massVarName = mName; }
  inline void  setMassVarSet(   bool set = true  ) { _massVarSet = set; }

  inline void addMassVar(float mmin = -1, float mmax = -1 ) {
    if( mmin > 0 ) _fitMassMin = mmin;
    if( mmax > 0 ) _fitMassMax = mmax;   
    RooRealVar *mass = new RooRealVar(massVarName().c_str(),"mass gg",_fitMassMin,_fitMassMax);
    _variables.add(*mass);
    _massVarSet = true;
  }


  //-- increase lumi (only sig or sig+bkg) for prospects
  void scaleLumiExp( float scaleSig, float scaleBkg = -1);

private:
  bool _massVarSet;
  bool _addSyst;
  TCut _blindingCut;
  TCut _mainCut;
  std::vector<TCut> _categories;
  std::vector<TString> _categoriesNames;
  std::vector<int>  _polynomialOrder;
  std::vector<std::vector<RooDataSet*> > _datasetSignalWoXS;
  std::vector<std::vector<RooDataSet*> > _datasetSignalWiXS;
  std::vector<std::vector<double> > _signalSigmaEff;
  std::vector<RooDataSet*> _dataset;

  int _bkgModel;
  std::string _massVarName;
  std::string _minitreeName;
  std::string _hlFactoryName;
  std::string _plotDirectory;
  RooStats::HLFactory    *_hlf;
  RooArgSet _variables;
  
  float _fitMassMin;
  float _fitMassMax;

  float _altSigAccCorr;

  /// unblinding flag
  bool _unblind;
  std::vector<bool>   _sigFitDone;
  std::vector<std::string> _pdfSigName;
  std::string         _pdfBkgName;

  /// simultaneous fit of several categories (for check)
  std::vector<float>            _nBkg;
  std::vector<float>            _nBkgExp;
  std::vector<float>            _nBkgPol;
  std::vector<float>            _slopeBkg;
  std::vector<float>            _texpBkg;
  std::vector<std::vector<float> >   _polBkg;
  std::vector<RooRealVar*>      _nSigExp;
  std::vector<RooSimultaneous*> _simPdf;
  std::vector<std::vector<RooAbsPdf*> > _bkgStudyPdf;
  RooCategory *_category;
  double simFitPreparePdf(int sHyp);

  void preparePseudoData(int sHyp, std::vector<RooDataSet*> &pseudoData_c);
  RooDataSet * prepareCombinedData( std::vector<RooDataSet*> &data );
  void preparePseudoDataBinned(int sHyp, std::vector<RooDataHist*> &pseudoData_c);
  RooDataHist * prepareCombinedData( std::vector<RooDataHist*> &data );

  void resetBackgroundPdf(void);
  float _scaleSig, _scaleBkg;
  
};

inline void MiniTreeFitter1D::setPlotDirectory( std::string dir ) { 
  _plotDirectory = dir + "/";
  gSystem->Exec( std::string("mkdir -p " + _plotDirectory ).c_str() );
}

MiniTreeFitter1D::MiniTreeFitter1D( std::string hlFactory ) {
  _massVarSet  = false;
  _massVarName = "mass";

  _addSyst = false;
  _altSigAccCorr = 1;

  _hlf = new RooStats::HLFactory("HLFactory", hlFactory.c_str(), false);
  _minitreeName = "HToGG";
  _unblind = false;

  _blindingCut = "mass < 115 || mass > 150";

  _fitMassMin = 100;
  _fitMassMax = 180;

  _scaleSig = 1;
  _scaleBkg = 1;

  /// add 2 basic variables always required
  //  RooRealVar *scaleXS = new RooRealVar("scaleXS","scale",0,10000);
  RooRealVar *wAccp = new RooRealVar("wei"   ,"acceptance weight "  ,0,100);
  RooRealVar *wXsec = new RooRealVar("wei_xs","cross section weight",0,100);
  _variables.add(*wAccp);
  _variables.add(*wXsec);

  _category = 0;
  _bkgModel = -1;
}

MiniTreeFitter1D::~MiniTreeFitter1D(void) {
  delete _hlf;

  for( unsigned i = 0 ; i < _dataset.size(); ++i ) delete _dataset[i];

  for( unsigned i = 0 ; i < _datasetSignalWoXS.size(); ++i )
  for( unsigned j = 0 ; j < _datasetSignalWoXS[i].size(); ++j )
    delete _datasetSignalWoXS[i][j];

  for( unsigned i = 0 ; i < _datasetSignalWiXS.size(); ++i )
  for( unsigned j = 0 ; j < _datasetSignalWiXS[i].size(); ++j )
    delete _datasetSignalWiXS[i][j];
 

  for( unsigned i = 0 ; i < _simPdf.size() ; ++i ) delete _simPdf[i];
  for( unsigned i = 0 ; i < _nSigExp.size(); ++i ) delete _nSigExp[i];

  delete _category;
 
  
}

#endif
