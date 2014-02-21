#include "Analysis/src/ToyAnalysis.hh"
 
#include <cassert>
#include <algorithm>
#include <sstream>
#include <fstream>
using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/CutUtils.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/selectors/ElectronSelector.hh"
#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/LeptonSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"

ClassImp( ToyAnalysis )

typedef vector<float> vectorFloat; 


ToyAnalysis::ToyAnalysis( Sample& sample, EventManager& manager ) : 
  SampleAnalysis( "Toy", sample, manager )
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---  Toy preselection analysis     ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _hlt = false; // apply HLT selection
  _nDebug = 0;
  hist_mll = new TH1F("m_ll", "m_ll", 100, 0, 200);
}

ToyAnalysis::~ToyAnalysis()
{
}

void
ToyAnalysis::bookHistograms()
{
//  defineTemplate( "Eproj", 400, -200, 200 );
//  defineTemplate( "Eproj2D", 200, -200, 200, 200, -200, 200 );
//  defineTemplate( "METVsMT", 80, 0, 80, 80, 0, 160 );
//  defineTemplate( "sigMET",  200, 0, 20    );
//  defineTemplate( "balance", 200, 0, 20    );
//  defineTemplate( "jmult", 20, -0.5, 19.5 );
//  defineTemplate( "LP", 100, -2., 2. );
  defineTemplate("m_ll", 100, 0, 200);
  tm.setTree("testNTuple", "");
  
}

void
ToyAnalysis::writeHistograms()
{
//  for( map< string, TH1* >::const_iterator it=h_.begin(); 
//       it!=h_.end(); it++ )
//    {
//      it->second->Write();
//    }
//  for( map< string, TH2* >::const_iterator it=h2_.begin(); 
//       it!=h2_.end(); it++ )
//    {
//      it->second->Write();
//    }  
//
  hist_mll->Write();
  tm.flush();
}

bool
ToyAnalysis::analyzeEvent()
{ 
  
  // build lists of leptons and prepare the MC matching map
  buildLeptonLists();


  const CandList& Electrons = _e.electrons();
  const CandList& Muons = _e.muons();
  Candidate* ZCand;
  int nElectrons, nMuons;
  nElectrons = Electrons.size();
  nMuons = Muons.size();
  double m_ll;
  if ( nElectrons >= 2 ) {
	  // create Z candidate with 2 electrons
	  ZCand = Candidate::create(Electrons[0], Electrons[1]);
	  cout << "ZCand->mass(): " << ZCand->mass() << endl;
	  m_ll = ZCand->mass();
	  mll_tree = m_ll;
	  tm.add<double>("m_ll", &mll_tree);
	  fill("m_ll", "ee", m_ll, "");

  }
  else if ( nMuons >= 2 ) {
	  // create Z candidate with 2 muons
	  ZCand = Candidate::create(Muons[0], Muons[1]);
	  cout << "ZCand->mass(): " << ZCand->mass() << endl;
	  m_ll = ZCand->mass();
	  mll_tree = m_ll;
	  tm.add<double>("m_ll", &mll_tree);
	  fill("m_ll", "mumu", m_ll, "");
  }
  else {
	  cout << "no Z candidate in this event" << endl;
	  m_ll = -1;
	  return false;
  }
  hist_mll->Fill(m_ll);

  // Retrieve jets.
  const CandList& Jets = _e.jetList( EventManager::kPfJet);
  Candidate* jet0;
  if ( Jets.size()!=0 ) {
	 jet0 = Jets[0];
  }
  double jet0_pt;
  if (jet0) {
	jet0_pt = jet0->pt();
	if( jet0_pt > 30 ) {
		tm.add<double>("jet0_pt", &jet0_pt);
	}
  }
  //for(int unsigned ij=0;ij<Jets.size();ij++) { 
  //  Candidate* jet = Jets[ij];
  //  if((fabs(jet->eta())<1.4442 ||
  //  (fabs(jet->eta())<2.5 && fabs(jet->eta())>1.56 ) ) )
  //    {
  //  if(pttmp<jet->pt())
  //    { theJet_=jet; pttmp = jet->pt(); }
  //    }
  //}

  // Retrieve met candidates.
  //Candidate* met_  = _e.met( EventManager::kPfMet);
 

  // stat histograms
  //fillStatHistograms();

  // debug printouts
  //debugPrintouts( cout );


  //return nSelected>0;
  return true;
}

