// fastjet stuff
#include "DISCambridgePlugin.hh"
#include "fastjet/NNH.hh"
#include <cassert>

// strings and streams
#include <sstream>
#include <limits>

// for kt2
#include <cmath>
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

namespace contrib{

// --------------------------------------------
// class to help run a DIS Cambridge algorithm
// it contains technical parameters, i.e. the p() power
// dij = min(Ei^2p, Ej^2p)*(1-cij)
// and the kT2() cut above which ij recombinations are not performed
class DISCamInfo{
  public: 
    DISCamInfo(double p, double kT2, int pz_beam_sign){
      _p= p; _kT2=kT2; _pz_beam_sign= pz_beam_sign;
    }

    double p() {return _p;}
    double kT2() {return _kT2;}
    int pz_beam_sign() {return _pz_beam_sign;}

  private:
    double _p, _kT2; 
    int _pz_beam_sign;

};


class DISCamBriefJet {
  public:
    // Jet properties used to compute distances
    void init(const PseudoJet & jet, DISCamInfo * info) {
        double norm = 1.0/sqrt(jet.modp2());
        _nx = jet.px() * norm;
        _ny = jet.py() * norm;
        _nz = jet.pz() * norm;
        _p = info->p();
        // _E2p = _p != 0 ? 1 : std::pow(jet.E(), 2*_p);
        _E2p = (_p == 0) ? 1 : std::pow(jet.E(), 2*_p);
        _pz_beam_sign = info->pz_beam_sign();
    }

    // Distance between *this and another *jet
    double distance(const DISCamBriefJet * jet) const{
        double dij = 1 - _nx*jet->_nx
                       - _ny*jet->_ny
                       - _nz*jet->_nz;
        dij = dij * min(_E2p, jet->_E2p);              
        return dij;
    }

    // Distance between *this and the beam (notice the beam is assumed to have pz>0)
    double beam_distance() const {
	    return (1 - _pz_beam_sign * _nz)*_E2p;
    }

    private:
      // normalised vector of three-momentum 
      double _nx, _ny, _nz;
      // Power of the energy to be used in the algorithm
      // and E^{2p}
      double _p, _E2p;

      int _pz_beam_sign; 
  };

//-----------------------------------------------------------------
  std::string DISCambridgePlugin::description () const {
    std::ostringstream desc;
    desc << "DISCambridge with kT2 cut= " << kT2() 
         << " and p = " << p() 
         << " and pz_beam_sign = " << pz_beam_sign(); 
    assert(p()>=0 && "p should be larger or equal than 0");
    return desc.str();
  }

//-----------------------------------------------------------------
// Clustering routine
  void DISCambridgePlugin::run_clustering(fastjet::ClusterSequence & cs) const{
    int njets = cs.jets().size();
     
    DISCamInfo info(p(), kT2(), pz_beam_sign());
    NNH<DISCamBriefJet,DISCamInfo> nnh(cs.jets(), &info);

    while (njets > 0) {

      // std::cout << "Begin with " << njets << " input jets with momentum " << std::endl;
      // for(int l=0; l< cs.jets().size(); l++){
        // std::cout<<"p"<<l<<"= "<< cs.jets()[l].px()<<" "<<cs.jets()[l].py()<<" "<<cs.jets()[l].pz()<<"\n";
      // }

      int i, j, k;
      double dij = nnh.dij_min(i, j);

      // std::cout<<"Smallest distance is between i and j = "<<i<<" "<<j;
      // std::cout<<" with dij = "<<vij<<"\n";

      // double vijrec;
      // if(j<0){
      //   vijrec = 1-cs.jets()[i].pz()/cs.jets()[i].modp();
      // }else{
      //   vijrec = 1-(cs.jets()[i].pz()*cs.jets()[j].pz()+ cs.jets()[i].py()*cs.jets()[j].py()+cs.jets()[i].px()*cs.jets()[j].px())/(cs.jets()[i].modp()*cs.jets()[j].modp());
      // }
      // std::cout<<" Our recomputed value is "<<vijrec<<"\n";
    
      // double diB = 2 * cs.jets()[i].beam_distance() * cs.jets()[i].E() * cs.jets()[i].E();  
      // double djB = 2 * cs.jets()[j].beam_distance() * cs.jets()[j].E() * cs.jets()[j].E(); 


      // Branchings above the kT2() cut are not clustered
      // Obtain the vij (= kT measure) from dij
      // dij = min(Ei^2p, Ej^2p) (1-cij)
      // vij = 2 min(Ei^2, Ej^2) (1-cij) = 2 min(Ei, Ej)^2 (1-cij) 
      // Since p>0 min(Ei^2p, Ej^2p) = min(Ei, Ej)^2p, thus
      // vij = 2 min(Ei, Ej)^{2-2p +2p} (1-cij) = 2 min(Ei, Ej)^{2-2p} dij
      double Emin = (j>=0) ? std::min(cs.jets()[i].E(), cs.jets()[j].E()) : cs.jets()[i].E();
      double vij = 2 * dij * pow(Emin, 2 - 2 * p()); //<-- kT^2 of the branching, to be compared against the kT2() cut

      // Veto an ij recombination if above the kT2() cut, and instead force beam recombination for the softest pron (i.e. tag it as jet)
      if (j>=0) {

        if (vij > kT2()){
          // std::cout << "VETO: swapping i and j, setting j=-1" << std::endl;
	        if (cs.jets()[i].E() > cs.jets()[j].E()) std::swap(i,j);
            j = -1;
        } 
      } 

      if (j >= 0) {
        cs.plugin_record_ij_recombination(i, j, dij, k);
        nnh.merge_jets(i, j, cs.jets()[k], k);
      } else {
        cs.plugin_record_iB_recombination(i, dij);
        nnh.remove_jet(i);
      }
      njets--;
      
    }
  }
}
  FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
