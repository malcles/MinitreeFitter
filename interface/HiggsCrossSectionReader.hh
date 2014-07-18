#ifndef higgsXSreader_HH_
#define higgsXSreader_HH_

#include <string>
#include <map>

std::map<float,float> ReadHiggsCrossSecionFile(    std::string file );
std::map<float,float> ReadHiggsBranchingRatioFile( std::string file );
class SMHiggsCrossSection {
 public:
  SMHiggsCrossSection( void );
  ~SMHiggsCrossSection(void);

  void is7TeV( void ) { _xsec_at_nTeV =  7; }
  void is8TeV( void ) { _xsec_at_nTeV =  8; }
  void is14TeV(void ) { _xsec_at_nTeV = 14; }

  float HiggsSMxsec_ggh( float mass );
  float HiggsSMxsec_vbf( float mass );
  float HiggsSMxsec_wh(  float mass );
  float HiggsSMxsec_zh(  float mass );
  float HiggsSMxsec_tth( float mass );
  float HiggsBR(  float mass );

 private:
  int     _xsec_at_nTeV;
  std::string _file_7TeV_ggh;
  std::string _file_8TeV_ggh;
  std::string _file_14TeV_ggh;
  std::string _file_7TeV_vbf;
  std::string _file_8TeV_vbf;
  std::string _file_14TeV_vbf;
  std::string _file_7TeV_wh;
  std::string _file_8TeV_wh;
  std::string _file_14TeV_wh;
  std::string _file_7TeV_zh;
  std::string _file_8TeV_zh;
  std::string _file_14TeV_zh;
  std::string _file_7TeV_tth;
  std::string _file_8TeV_tth;
  std::string _file_14TeV_tth;

  std::string _file_br;

  std::map<float,float> _xsec_ggh;
  std::map<float,float> _xsec_vbf;
  std::map<float,float> _xsec_wh;
  std::map<float,float> _xsec_zh;
  std::map<float,float> _xsec_tth;
  std::map<float,float> _br;


};


#endif
