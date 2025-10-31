#ifndef __FASTJET_CONTRIB_DISCAMBRIDGEALGORITHM_HH__
#define __FASTJET_CONTRIB_DISCAMBRIDGEALGORITHM_HH__


#include <fastjet/internal/base.hh>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{
  
    class DISCambridgePlugin : public JetDefinition::Plugin {
    public:

        // Constructor for the DIS Cambridge Plugin class - taking parameters kT2 and p
	DISCambridgePlugin(double kT2, double p, int pz_beam_sign): _kT2(kT2), _p(p), _pz_beam_sign(pz_beam_sign) {}

                /// Copy constructor (member-wise copy to avoid relying on implicitly-declared assignment)
                DISCambridgePlugin (const DISCambridgePlugin & plugin)
                        : _kT2(plugin._kT2), _p(plugin._p), _pz_beam_sign(plugin._pz_beam_sign) {}

        // the things that are required by base class
        virtual std::string description () const;
        virtual void run_clustering(ClusterSequence &) const;

        // kT2() cut
        double kT2() const {return _kT2;}
        
	// Power of energy scale
        double p() const{return _p;}

	// Beam direction
        int pz_beam_sign() const{return _pz_beam_sign;}

	/// Returns a dummy jet radius as not required for this algorithm
        virtual double R() const {return 1.0;}

        /// Avoids warning when user requests "exclusive" jets from cluster sequence
        virtual bool exclusive_sequence_meaningful() const {return true;}

        /// returns false because this plugin is not intended for spherical
        /// geometries (i.e. it's a DIS algorithm).
        virtual bool is_spherical() const override {return false;}

    private:
	double _kT2;
        double _p;
        int _pz_beam_sign;
    };

} //namespace contrib


FASTJET_END_NAMESPACE

#endif 
