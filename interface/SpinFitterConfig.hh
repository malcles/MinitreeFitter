#ifndef spinFitter__hh__
#define spinFitter__hh__

#include <TCut.h>

#include <iostream>
#include <vector>
#include <string>

class SpinFitterConfig {
public:
  SpinFitterConfig() {}
  SpinFitterConfig(std::string configname);
  ~SpinFitterConfig(void) {}
  void setConfigFile(std::string configname) { 
    _config = configname; setup();
  }
  std::vector<float> getCosThetaStarBoundaries(void)                const {return _cosThetaStarBoundaries; }
  std::vector<float> getDiphoDiscriBoundaries(unsigned cosThetaBin) const {
    return cosThetaBin<_cosDiphoDiscriBoundaries.size()? 
      _cosDiphoDiscriBoundaries[cosThetaBin]:std::vector<float>(); 
  }
  std::vector<TCut>          getCosThetaCuts(void) const;
  std::vector<std::vector<TCut> > getDiphoDiscriCuts(void) const;
  std::vector<TCut>          getDiphoDiscriCuts(unsigned cosThetaBin) const { 
    return cosThetaBin<getDiphoDiscriCuts().size()?
      getDiphoDiscriCuts()[cosThetaBin]:std::vector<TCut>(); 
  }

  float AccTot(std::string rootfile, std::string treename = "HToGG") const;
  
private: 
  std::string _config;
  std::vector<float>               _cosThetaStarBoundaries;
  std::vector<std::vector<float> > _cosDiphoDiscriBoundaries;
  void setup();
};


#endif
