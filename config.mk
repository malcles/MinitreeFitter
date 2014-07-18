LIBRARYNAME := MinitreeFitter


ifeq "$(shell uname -n | head -c6)" "lxplus"
  BOOST = /afs/cern.ch/cms/slc5_amd64_gcc434/external/boost/1.47.0
  BOOST_INCLUDES := 
  BOOST_LIBS :=  -lboost_program_options
else
  BOOST_INCLUDES := 
  BOOST_LIBS := -lboost_program_options
endif

BOOST_LIBS += -lRooFit -lRooFitCore -lRooStats