#ifndef ToyAnalysis_h
#define ToyAnalysis_h

#include <vector>
#include <map>
#include <string>

class TH1;
class TH2;

#include "Analysis/core/SampleAnalysis.hh"

#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/ElectronSelector.hh"

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
  //XXX TH1* hist_mll;
  //XXX double mll_tree;
  float getRapidity(const Candidate* zCandidate);

  MuonSelector MuonSel;
  ElectronSelector ElectronSel;

  ClassDef( ToyAnalysis, 0 )  
};

#endif

     
