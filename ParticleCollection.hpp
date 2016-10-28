#ifndef _NANOBD_PARTICLECOLLECTION_HPP
#define _NANOBD_PARTICLECOLLECTION_HPP

#include <vector>
#include "Particle.hpp"
#include "BondList.hpp"
#include "PairList.hpp"
#include "CellCollection.hpp"
#include "config.hpp"
#include <gsl/gsl_matrix.h>

class ParticleCollection : public std::vector<Particle> {
  public:
    ParticleCollection() : std::vector<Particle>() {};
   
    void initialize(double& sidex_, double& sidey_, double& sidez_, const double a0_);
   
    void updateBondList(const PairList pairs_, double sidex_,  double sidey_, double sidez_, double dt_,
                        BondList& bonds_, const BondList unbreakable_,
                        double xcp_[], double Lcp_[], long long t_,
                        VectorLL& life_times_, size_t num_bond_types_);
    void computeExcludedVolumeForces(const PairList pairs_, double sidex_, double sidey_, 
                                     double sidez_, bool remove_overlap_);

    void computeForcesAndTorques(const PairList pairs_,double sidex_,  double sidey_, double sidez_, 
                                size_t max_site_type);
//    void computeNonPatchForcesAndTorques(const PairList pairs_, double sidex_,  double sidey_, double sidez_,size_t max_site_type);
    void applyRepulsiveForce(double sidex_,  double sidey_, double sidez_);
//    void applyRepulsiveForce2(double sidex_,  double sidey_, double sidez_);
//    void computeAngles(long long t_, double dt_, const PairList pairs_, double sidex_,  double sidey_, double sidez_ ,size_t max_site_type);
    void updateDiffusionTensor(const PairList pairs_,double sidex_,  double sidey_, double sidez_ ,
                               const double a0_, gsl_matrix* Dij_);
    void computeSpringForcesAndTorques(double sidex_,  double sidey_, double sidez_, BondList bonds, double Lcp_[], double k_spring_);
    void updatePositionsAndOrientations(double sidex_,  double sidey_, double sidez_, double dt_, int verlet_, double rc_,
                                        double rs_, bool& skin_violation_, const double a0_,
                                        gsl_matrix* Dij_, bool hydrodynamics_);
    void printPositions(long long t_) const;
    void printSiteCollection(long long t_) const;
    void printMSD(long long t_, double sidex_,  double sidey_, double sidez__, double dt_) const;
    void printCosTheta(long long t_, double dt_) const;
    void printBonds(long long t_, double dt_, const BondList& bonds_) const;
    void writeBondStats(long long t_, double dt_, const BondList bonds_);
//    void writeContactsAB(long long t_, double dt_, double sidex_,  double sidey_, double sidez_);
 //   void writeContactsAAorBB(long long t_, double dt_,double sidex_,  double sidey_, double sidez_);
  //  void writeContactsALL(long long t_, double dt_,double sidex_,  double sidey_, double sidez_);
    void writeClusterStats(long long t_, double dt_, const BondList bonds_) const;
    void writeClusterSizes(long long t_, double dt_, const BondList bonds_, bool composition_) const;
    void writeBondLifeTimes(VectorLL& life_times_, double dt_) const;
    void buildPairs(PairList& pairs_, CellCollection& cells_,const double sidex_,  const double sidey_, const double sidez_ , 
                     const int verlet_, double& rs_, const double rc_, bool& skin_violation_, const double rs_next) const;
    void readScript(const char* s_, long long& step_, long long& steps_, long long& freqWrite_, double& dt_,
                    double& sidex_, double& sidey_, double& sidez_, double& k_spring_, int& verlet_, double& rc_,
                    double& rs_, BondList& bonds_, unsigned int& seed_, bool& restart_,
                    BondList& unbreakable_, double xcp_[], double Lcp_[],
                    size_t& max_site_type_, const double a0_);
    void readConfiguration(const char* s_, long long& step_, double& sidex_, double& sidey_, double& sidez_, BondList& bonds_,
                           BondList& unbreakable_, const double a0_);
    void writeConfiguration(long long t_, long long steps_, long long freqWrite_, double dt_, double sidex_, double sidey_, double sidez_, double k_spring_, BondList bonds_, BondList unbreakable_,
                            const double a0_) const;
};

#endif
