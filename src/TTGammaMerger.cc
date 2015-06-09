// -*- C++ -*-
//
// Package:    TTGammaMerger
// Class:      TTGammaMerger
// 
/**\class TTGammaMerger TTGammaMerger.cc MyPackage/TTGammaMerger/src/TTGammaMerger.cc

 Description: Removes Signal Overlap in TTbar sample

 Implementation:
// Sort out events, that have been simulated with ttgamma matrix element.

*/
//
// Original Author:  Heiner Tholen
//         Created:  Wed May 23 20:38:31 CEST 2012
// $Id: TTGammaMerger.cc,v 1.7 2013/06/27 14:30:17 htholen Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <TH1D.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class TTGammaMerger : public edm::EDFilter {
   public:
      explicit TTGammaMerger(const edm::ParameterSet&);
      ~TTGammaMerger();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void printParticle(const reco::GenParticle* p, std::ostringstream &out);
      virtual void printParticles(std::vector<const reco::GenParticle*> &v, std::ostringstream &out);
      virtual void findMothers(const reco::GenParticle* p, std::vector<const reco::GenParticle*> &moms);
      virtual void findDaughters(const reco::GenParticle* p,std::vector<const reco::GenParticle*> &all);
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      const double ptCut_;
      const double drCut_;
      const double etaCut_;
      TH1D *etKickedPhotons_;
      TH1D *etSurvivingPhotons_;
      TH1D *etAllPhotons_;
      TH1D *etaKickedPhotons_;
      TH1D *etaSurvivingPhotons_;
      TH1D *etaAllPhotons_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TTGammaMerger::TTGammaMerger(const edm::ParameterSet& iConfig) :
    ptCut_(iConfig.getParameter<double>("ptCut")),
    drCut_(iConfig.getParameter<double>("drCut")),
    etaCut_(iConfig.getParameter<double>("etaCut"))
{
    edm::Service<TFileService> fs;
    etaKickedPhotons_     = fs->make<TH1D>("etaKickedPhotons",    ";photon E_{T} / GeV;number of photons", 80, -4., 4.);
    etaSurvivingPhotons_  = fs->make<TH1D>("etaSurvivingPhotons", ";photon E_{T} / GeV;number of photons", 80, -4., 4.);
    etaAllPhotons_        = fs->make<TH1D>("etaAllPhotons",       ";photon E_{T} / GeV;number of photons", 80, -4., 4.);
    etKickedPhotons_      = fs->make<TH1D>("etKickedPhotons",     ";photon E_{T} / GeV;number of photons", 70, 0., 700.);
    etSurvivingPhotons_   = fs->make<TH1D>("etSurvivingPhotons",  ";photon E_{T} / GeV;number of photons", 70, 0., 700.);
    etAllPhotons_         = fs->make<TH1D>("etAllPhotons",        ";photon E_{T} / GeV;number of photons", 70, 0., 700.);

    produces<std::vector<reco::GenParticle> >("signalPhotons");
}


TTGammaMerger::~TTGammaMerger()
{
}


//
// member functions
//

void
TTGammaMerger::printParticle(const reco::GenParticle* p, std::ostringstream &out)
{
    char buf[256];
    snprintf(buf, 256,
           " %5d | %2d | %7.3f %10.3f %6.3f | %10.3f %10.3f %10.3f %8.3f |\n",
           p->pdgId(),
           p->status(),
           p->pt(),
           p->eta(),
           p->phi(),
           p->px(),
           p->py(),
           p->pz(),
           p->mass()
          );
    out << buf;
}

void
TTGammaMerger::printParticles(std::vector<const reco::GenParticle*> &v, std::ostringstream &out)
{
    for (unsigned i = 0; i < v.size(); ++i) printParticle(v.at(i), out);
}

void
TTGammaMerger::findMothers(
    const reco::GenParticle* p,
    std::vector<const reco::GenParticle*> &moms
) {
    // only partons are interesting:
    if (abs(p->pdgId()) > 21) return;

    // already in vector of moms?
    bool new_particle = true;
    for (unsigned i = 0; i < moms.size(); ++i) {
        if (moms[i] == p) {
            new_particle = false;
            break;
        }
    }
    if (new_particle) {
        moms.push_back(p);

        // go recursive
        for (unsigned i = 0; i < p->numberOfMothers(); ++i) {
            findMothers((const reco::GenParticle*) p->mother(i), moms);
        }
    }
}

void
TTGammaMerger::findDaughters(
    const reco::GenParticle* p,
    std::vector<const reco::GenParticle*> &all
) {
    // walk over daughters
    for (unsigned i = 0; i < p->numberOfDaughters(); ++i) {
        const reco::GenParticle* da = (const reco::GenParticle*) p->daughter(i);
        int abs_pdg = abs(da->pdgId());

        // take only tops, bs and Ws
        if (abs_pdg == 6 || abs_pdg == 24) {
            all.push_back(da);
            findDaughters(da, all);
        }
  else if (abs_pdg < 6 ||(abs_pdg >10 && abs_pdg < 17)){
      all.push_back(da);
  }
    }
}

// ------------ method called on each new Event  ------------
bool
TTGammaMerger::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using reco::GenParticle;
    using reco::deltaR;

    ostringstream out("TTGammaMerger");
    out << "    ID |Stat|    pt       eta     phi   |     px         py         pz        m     |\n";

    Handle<vector<GenParticle> > gens;
    iEvent.getByLabel(InputTag("genParticles"), gens);

    std::vector<reco::GenParticle>* signalPhotons = new std::vector<reco::GenParticle>();
    std::auto_ptr<std::vector<reco::GenParticle> > pOut(signalPhotons);

    /////////////////////////////////////////////////// find core particles ///
    vector<const GenParticle*> all;

    // find first top quarks as a starting point
    GenParticle local_top;
    GenParticle local_topBar;
    const GenParticle* top      = 0;
    const GenParticle* topBar   = 0;
    for (vector<GenParticle>::const_iterator i = gens->begin(); i != gens->end(); ++i){
         if (!top && i->pdgId() ==  6) {
            local_top = *i;
            top = &local_top;
         }
         if (!topBar && i->pdgId() == -6) {
            local_topBar = *i;
            topBar = &local_topBar;
         }
         if (top && topBar) break;
    }
    out << "top, topBar\n";
    printParticle(top, out);
    printParticle(topBar, out);

    // find initial state particles
    findMothers(top, all);
    findMothers(topBar, all);
    out << "top, topBar and moms\n";
    printParticles(all, out);

    // find relevant final states
    findDaughters(top, all);
    findDaughters(topBar, all);
    out << "top, topBar, moms and daughters\n";
    printParticles(all, out);

    ///////////////////////////////////////////////// find relevant photons ///
    vector<const GenParticle*> photons;
    for (unsigned i = 0; i < all.size(); ++i) {
        for (unsigned j = 0; j < all.at(i)->numberOfDaughters(); ++j) {
            const GenParticle* daughter = (const GenParticle*) all.at(i)->daughter(j);
            if (abs(daughter->pdgId()) == 22) {
                 photons.push_back(daughter);
                 etAllPhotons_->Fill(daughter->et());
                 etaAllPhotons_->Fill(daughter->eta());
            }
        }
    }
    out << "photons\n";
    printParticles(photons, out);

    // put signal photons to event
    for (unsigned i = 0; i < photons.size(); ++i) signalPhotons->push_back(*photons.at(i));
    iEvent.put(pOut, "signalPhotons");

    /////////////////////////////////////////////////////// find b's, lightquarks and leptons ///
    const GenParticle* b    = 0;
    const GenParticle* bbar = 0;
    
    vector<const GenParticle*> lightquarks;
    vector<const GenParticle*> leptons;
    vector<const GenParticle*> neutrinos;
    vector<const GenParticle*> quarks;
    vector<const GenParticle*> legs;

    // b's are used to check for correct decay
    for (int i = all.size() - 1; i > -1; --i) {
        const GenParticle* p = all.at(i);
        //int stat = p->status();
        int pdg  = p->pdgId();
  int abs_pdg = abs(pdg);
        if (!b    && pdg ==  5) b    = p;
        if (!bbar && pdg == -5) bbar = p;
        if (b && bbar) break;
        if (abs_pdg < 5) lightquarks.push_back(p);
        else if (abs_pdg == 11 || abs_pdg == 13 || abs_pdg == 15 ) leptons.push_back(p);
    }

    if (!(b && bbar)) return true; // top did not decay to W+b: return

    vector<const GenParticle*> bquarks;
    bquarks.push_back(b);
    bquarks.push_back(bbar);
    out << "bquarks (b, bbar)\n";
    printParticles(bquarks, out);

    //////////////////Check for lepton, b and jet cuts //////////////
    quarks.insert(quarks.end(),bquarks.begin(),bquarks.end());
    quarks.insert(quarks.end(),lightquarks.begin(),lightquarks.end()); 
   
    for (unsigned i = 0; i < quarks.size();++i){
  const GenParticle* quark = quarks.at(i);
  if (quark->pt() > 20) return true;
        for (unsigned j = 0; j < quarks.size();++j ){
      const GenParticle* leg = quarks.at(j);
      if (deltaR(*quark, *leg) < 0.5) return true;
  }
  for (unsigned k= 0; k<leptons.size();k++){
      const GenParticle* lepton = leptons.at(k);
      if (deltaR(*quark, *lepton) < 0.5) return true;
      }
    }

    /////////////////////////////// sort out fails (must fulfill all three cuts) ///
    //////////////////////Only Works if All Photon delta R Cuts are the same ///////
    
    legs.insert(legs.end(),quarks.begin(),quarks.end());
    legs.insert(legs.end(),leptons.begin(),leptons.end());
    bool foundNoSignalPhoton = true;
    for (unsigned i = 0; i < photons.size(); ++i) {
        const GenParticle* photon = photons.at(i);
        if (photon->pt() > ptCut_ && photon->eta() < etaCut_) {
            bool closeToLeg = false;
            for (unsigned j = 0; j < legs.size(); ++j) {
                const GenParticle* leg = legs.at(j);
                if (deltaR(*photon, *leg) < drCut_) {
                    closeToLeg = true;
                }
            }
            if (!closeToLeg) {
                foundNoSignalPhoton = false;
                out << "<TTGammaMerger>: removing Event! "
                    << "Photon pt < ptCut: (" << photon->pt() << " > " << ptCut_
                    << ") and no deltaR to a leg smaller than " << drCut_ << endl;
                break;
                
            }
        }
    }

    if (foundNoSignalPhoton) {
        for (unsigned i = 0; i < photons.size(); ++i) {
            etSurvivingPhotons_->Fill(photons.at(i)->et());
            etaSurvivingPhotons_->Fill(photons.at(i)->eta());
        }
        LogInfo("TTGammaMerger") << out.str();
     } else {
        for (unsigned i = 0; i < photons.size(); ++i) {
            etKickedPhotons_->Fill(photons.at(i)->et());
            etaKickedPhotons_->Fill(photons.at(i)->eta());
        }
        LogWarning("TTGammaMerger") << out.str();
    }
    return foundNoSignalPhoton;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTGammaMerger::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTGammaMerger::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
TTGammaMerger::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
TTGammaMerger::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
TTGammaMerger::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
TTGammaMerger::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTGammaMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TTGammaMerger);