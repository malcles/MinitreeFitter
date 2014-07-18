#ifndef setupreader_HH_
#define setupreader_HH_

#include <string>

class SetupReader {
 public:
  SetupReader( std::string setupfile );
  ~SetupReader(void);
  
  void removeSpaces(std::string &st) const;
  
  int    getInt(   std::string vname, bool &found ) const;
  float  getFloat( std::string vname, bool &found ) const;
  std::string getString(std::string vname, bool &found ) const;
  
  std::string setupFile(void) const {return _setupfile;}
  
 private:
  std::string _setupfile;
};

#endif
