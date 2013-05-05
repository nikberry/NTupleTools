#include "BristolAnalysis/NTupleTools/interface/BristolNTuple_Photons.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
 
BristolNTuple_Photons::BristolNTuple_Photons(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    vertexTag(iConfig.getParameter<edm::InputTag>("VertexTag")),
    particleFlowTag(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
    produces<std::vector<double> > (prefix + "Px" + suffix);
    produces<std::vector<double> > (prefix + "Py" + suffix);
    produces<std::vector<double> > (prefix + "Pz" + suffix);
    produces<std::vector<double> > (prefix + "Energy" + suffix);
    produces<std::vector<double> > (prefix + "EcalIso" + suffix);
    produces<std::vector<double> > (prefix + "HcalIso" + suffix);
    produces<std::vector<double> > (prefix + "HoE" + suffix);
    produces<std::vector<double> > (prefix + "TrkIso" + suffix);
    produces<std::vector<double> > (prefix + "SigmaIEtaIEta" + suffix);
    produces<std::vector<bool> > (prefix + "TrkVeto" + suffix);
    produces<std::vector<double> > (prefix + "SCseedEnergy" + suffix);
    produces<std::vector<double> > (prefix + "SCenergy" + suffix);
    produces<std::vector<double> > (prefix + "SCeta" + suffix);
    produces<std::vector<double> > (prefix + "SCphi" + suffix);
    produces<std::vector<double> > (prefix + "E3x3" + suffix);
    produces<std::vector<double> > (prefix + "E5x5" + suffix);
    produces<std::vector<double> > (prefix + "PfChargedIso03" + suffix);
    produces<std::vector<double> > (prefix + "PfPhotonIso03" + suffix);
    produces<std::vector<double> > (prefix + "PfNeutralIso03" + suffix);
    produces<std::vector<double> > (prefix + "HcalIso2012" + suffix);
    produces<std::vector<double> > (prefix + "HtowoE" + suffix);
    produces<std::vector<double> > (prefix + "ConvSafeEle" + suffix);
    
    isolator.initializePhotonIsolation(kTRUE);
    isolator.setConeSize(0.3);
       
}

void BristolNTuple_Photons::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::auto_ptr < std::vector<double> > px(new std::vector<double>());
    std::auto_ptr < std::vector<double> > py(new std::vector<double>());
    std::auto_ptr < std::vector<double> > pz(new std::vector<double>());
    std::auto_ptr < std::vector<double> > energy(new std::vector<double>());
    std::auto_ptr < std::vector<double> > ecalIso(new std::vector<double>());
    std::auto_ptr < std::vector<double> > hcalIso(new std::vector<double>());
    std::auto_ptr < std::vector<double> > hoe(new std::vector<double>());
    std::auto_ptr < std::vector<double> > trkIso(new std::vector<double>());
    std::auto_ptr < std::vector<double> > sigmaIetaIeta(new std::vector<double>());
    std::auto_ptr < std::vector<bool> > trkVeto(new std::vector<bool>());
    std::auto_ptr < std::vector<double> > SCseedEnergy(new std::vector<double>());
    std::auto_ptr < std::vector<double> > SCenergy(new std::vector<double>());
    std::auto_ptr < std::vector<double> > SCeta(new std::vector<double>());
    std::auto_ptr < std::vector<double> > SCphi(new std::vector<double>());
    std::auto_ptr < std::vector<double> > E3x3(new std::vector<double>());
    std::auto_ptr < std::vector<double> > E5x5(new std::vector<double>());
    std::auto_ptr < std::vector<double> > pfchargedIso(new std::vector<double>());
    std::auto_ptr < std::vector<double> > pfphotonIso(new std::vector<double>());
    std::auto_ptr < std::vector<double> > pfneutralIso(new std::vector<double>());
    std::auto_ptr < std::vector<double> > hcalIso2012(new std::vector<double>());
    std::auto_ptr < std::vector<double> > htowoe(new std::vector<double>());
    std::auto_ptr < std::vector<double> > convsafeele(new std::vector<double>());
    
    //-----------------------------------------------------------------

    edm::Handle < std::vector<reco::Photon> > photons;
    iEvent.getByLabel(inputTag, photons);
    
    //for pf reliso
    Handle<reco::VertexCollection>  vertexCollection;
    iEvent.getByLabel(vertexTag, vertexCollection);

    // All PF Candidate for alternate isolation
    Handle<reco::PFCandidateCollection> pfCandidatesH;
    iEvent.getByLabel(particleFlowTag, pfCandidatesH);
    const  PFCandidateCollection thePfColl = *(pfCandidatesH.product());

    unsigned int ivtx = 0;
    VertexRef myVtxRef(vertexCollection, ivtx);
   
   //for safe electron veto
   edm::Handle<reco::BeamSpot> bsHandle;
   iEvent.getByLabel("offlineBeamSpot", bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();

   edm::Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel("allConversions", hConversions);

   edm::Handle<reco::GsfElectronCollection> hElectrons;
   iEvent.getByLabel("gsfElectrons", hElectrons);
   
    if (photons.isValid()) {
        edm::LogInfo("BristolNTuple_PhotonsInfo") << "Total # Photons: " << photons->size();

        for (std::vector<reco::Photon>::const_iterator it = photons->begin(); it != photons->end(); ++it) {
            // exit from loop when you reach the required number of photons
            if (px->size() >= maxSize)
                break;
	  
	  
	    isolator.fGetIsolation(&*it, &thePfColl, myVtxRef,vertexCollection);
	    
	    //cout<<"PF  :  "<<isolator.getIsolationCharged()<<" : "<<isolator.getIsolationPhoton()<<" : "<<isolator.getIsolationNeutral()<<endl;
	    bool passelectronveto = !ConversionTools::hasMatchedPromptElectron(it->superCluster(), hElectrons, hConversions, beamspot.position());
	    //cout << "pass con: " << passelectronveto << endl;
	    
	    px->push_back(it->px());
            py->push_back(it->py());
            pz->push_back(it->pz());
            energy->push_back(it->energy());
            ecalIso->push_back(it->ecalRecHitSumEtConeDR04());
            hcalIso->push_back(it->hcalTowerSumEtConeDR04());
            hoe->push_back(it->hadronicOverEm());
            trkIso->push_back(it->trkSumPtHollowConeDR04());
            sigmaIetaIeta->push_back(it->sigmaIetaIeta());
            trkVeto->push_back(it->hasPixelSeed());
            SCseedEnergy->push_back(it->superCluster()->seed()->energy());
            SCenergy->push_back(it->superCluster()->energy());
            SCeta->push_back(it->superCluster()->eta());
            SCphi->push_back(it->superCluster()->phi());
            E3x3->push_back(it->e3x3());
            E5x5->push_back(it->e5x5());
	    pfchargedIso->push_back(isolator.getIsolationCharged());
            pfphotonIso->push_back(isolator.getIsolationPhoton());
	    pfneutralIso->push_back(isolator.getIsolationNeutral());
	    hcalIso2012->push_back(it->hcalTowerSumEtConeDR04() + (it->hadronicOverEm() - it->hadTowOverEm())*it->superCluster()->energy()/cosh(it->superCluster()->eta()));
	    htowoe->push_back(it->hadTowOverEm());
	    convsafeele->push_back(passelectronveto);
	}
    } else {
        edm::LogError("BristolNTuple_PhotonsError") << "Error! Can't get the product " << inputTag;
    }

    //-----------------------------------------------------------------
    // put vectors in the event
    iEvent.put(px, prefix + "Px" + suffix);
    iEvent.put(py, prefix + "Py" + suffix);
    iEvent.put(pz, prefix + "Pz" + suffix);
    iEvent.put(energy, prefix + "Energy" + suffix);
    iEvent.put(ecalIso, prefix + "EcalIso" + suffix);
    iEvent.put(hcalIso, prefix + "HcalIso" + suffix);
    iEvent.put(hoe, prefix + "HoE" + suffix);
    iEvent.put(trkIso, prefix + "TrkIso" + suffix);
    iEvent.put(sigmaIetaIeta, prefix + "SigmaIEtaIEta" + suffix);
    iEvent.put(trkVeto, prefix + "TrkVeto" + suffix);
    iEvent.put(SCseedEnergy, prefix + "SCseedEnergy" + suffix);
    iEvent.put(SCenergy, prefix + "SCenergy" + suffix);
    iEvent.put(SCeta, prefix + "SCeta" + suffix);
    iEvent.put(SCphi, prefix + "SCphi" + suffix);
    iEvent.put(E3x3, prefix + "E3x3" + suffix);
    iEvent.put(E5x5, prefix + "E5x5" + suffix);
    iEvent.put(pfchargedIso, prefix + "PfChargedIso03" + suffix);
    iEvent.put(pfphotonIso, prefix + "PfPhotonIso03" + suffix);
    iEvent.put(pfneutralIso, prefix + "PfNeutralIso03" + suffix);
    iEvent.put(hcalIso2012, prefix + "HcalIso2012" + suffix);
    iEvent.put(htowoe, prefix + "HtowoE" + suffix);
    iEvent.put(convsafeele, prefix + "ConvSafeEle" + suffix);
}
