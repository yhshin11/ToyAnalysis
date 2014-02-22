#include "Analysis/src/ToyAnalysis.hh"
 
#include <cassert>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <math.h>
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
SampleAnalysis( "Toy", sample, manager ),
MuonSel(MuonSelector::kPt20,MuonSelector::kTight),
ElectronSel(MuonSelector::kPt20,MuonSelector::kMedium)
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---  Toy preselection analysis     ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _hlt = false; // apply HLT selection
  _nDebug = 0;
  //XXX hist_mll = new TH1F("m_ll", "m_ll", 100, 0, 200);
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
  //XXX hist_mll->Write();
  tm.flush();
}

bool
ToyAnalysis::analyzeEvent()
{ 
  
  // build lists of leptons and prepare the MC matching map
  buildLeptonLists();


  const CandList& Electrons = _e.electrons();
  const CandList& Muons = _e.muons();
  int nElectrons, nMuons;
  nElectrons = Electrons.size();
  nMuons = Muons.size();
  Candidate *l1, *l2;
  // Lepton variables to store.
  bool boolTightIso_l1;
  int pdgCode_l1 = 0;
  double pt_l1, eta_l1, phi_l1, charge_l1;
  bool boolTightIso_l2;
  int pdgCode_l2 = 0;
  double pt_l2, eta_l2, phi_l2, charge_l2;


  Candidate* ZCand;
  // Z candidate variables to store.
  bool boolValidZ = false;
  double m_ll, pt_ll, phi_ll, y_ll;

  if ( nElectrons >= 2 ) {
	  l1 = Electrons[0];
	  l2 = Electrons[1];
	  // XXX fill("m_ll", "ee", m_ll, "");
	  boolValidZ = true;
	  boolTightIso_l1 = ElectronSel.accept(*l1);
	  boolTightIso_l2 = ElectronSel.accept(*l2);

  }
  else if ( nMuons >= 2 ) {
	  l1 = Muons[0];
	  l2 = Muons[1];
	  // XXX fill("m_ll", "mumu", m_ll, "");
	  boolValidZ = true;
	  boolTightIso_l1 = MuonSel.accept(*l1);
	  boolTightIso_l2 = MuonSel.accept(*l2);
  }
  else {
	  cout << "no Z candidate in this event" << endl;
	  boolValidZ = false;
  }
  if (boolValidZ) {
	pdgCode_l1 = l1->pdgCode();
	pt_l1 = l1->pt();
	eta_l1 = l1->eta();
	phi_l1 = l1->phi();
	charge_l1 = l1->charge();
	pdgCode_l2 = l2->pdgCode();
	pt_l2 = l2->pt();
	eta_l2 = l2->eta();
	phi_l2 = l2->phi();
	charge_l2 = l2->charge();
	// Create Z candidate with 2 leptons.
	ZCand = Candidate::create(l1, l2);
	m_ll = ZCand->mass();
	pt_ll = ZCand->pt();
	phi_ll = ZCand->phi();
	y_ll = getRapidity(ZCand);
  }
  tm.add<bool>("boolTightIso_l1", &boolTightIso_l1);
  tm.add<int>("pdgCode_l1", &pdgCode_l1);
  tm.add<double>("pt_l1", &pt_l1);
  tm.add<double>("eta_l1", &eta_l1);
  tm.add<double>("phi_l1", &phi_l1);
  tm.add<double>("charge_l1", &charge_l1);
  tm.add<bool>("boolTightIso_l2", &boolTightIso_l2);
  tm.add<int>("pdgCode_l2", &pdgCode_l2);
  tm.add<double>("pt_l2", &pt_l2);
  tm.add<double>("eta_l2", &eta_l2);
  tm.add<double>("phi_l2", &phi_l2);
  tm.add<double>("charge_l2", &charge_l2);

  tm.add<bool>("boolValidZ", &boolValidZ);
  tm.add<double>("m_ll", &m_ll);
  tm.add<double>("pt_ll", &pt_ll);
  tm.add<double>("phi_ll", &phi_ll);
  tm.add<double>("y_ll", &y_ll);
  //XXX hist_mll->Fill(m_ll);

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

double 
ToyAnalysis::getRapidity(const Candidate* zCandidate)
{
	if (zCandidate) {
		double pz = zCandidate->pz();
		double E = zCandidate->E();
		double rapidity = 0.5*log( (E + pz)/(E - pz) );
		return rapidity;
	}
	else return 0;
}
