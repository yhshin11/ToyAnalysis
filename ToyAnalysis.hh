#ifndef ToyAnalysis_h
#define ToyAnalysis_h

#include <vector>
#include <map>
#include <string>

class TH1;
class TH2;

#include "Analysis/core/SampleAnalysis.hh"

class ToyAnalysis : public SampleAnalysis
{
public:

  //  ToyAnalysis( Sample& sample, const std::string & collectionFileName );
  ToyAnalysis( Sample& sample, EventManager& manager );
  virtual ~ToyAnalysis();
 
public:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

private:
  // insert private members here
  TH1* hist_mll;
  double mll_tree;

  ClassDef( ToyAnalysis, 0 )  
};

#endif

     
