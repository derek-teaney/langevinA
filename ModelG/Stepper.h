#ifndef STEPPER_H
#define STEPPER_H

#include "ModelA.h"

class Stepper {
public:
  virtual bool step(const double &dt) = 0;
  virtual void finalize() = 0;
  virtual ~Stepper() = default;
};


/////////////////////////////////////////////////////////////////////////
//! Helper class to monitor the jumps in phi at each time step, and to
//! record the number of successful steps
class o4_stepper_monitor {
public:
  long int down_counts;    //!< number of steps downwards
  long int up_counts;      //!< number of upward trial steps
  long int up_yes;         //!< number of accepted upward steps
  long int up_no;          //!< number of rejected upward steps
  double down_dS;          //!< The change in action summed over downward steps
  double up_dS;            //!< The change in action summed over upward steps
  std::string description; //!< What this is monitoring

public:
  o4_stepper_monitor(const std::string &name = "monitor") : description(name) {
    down_counts = 0;
    up_counts = 0;
    up_yes = 0;
    up_no = 0;
    down_dS = 0;
    up_dS = 0;
  }
  void increment_up_yes(const double &deltaS) {
    up_counts++;
    up_yes++;
    up_dS += deltaS;
  }
  void increment_up_no(const double &deltaS) {
    up_counts++;
    up_no++;
    up_dS += 0.;
  }
  void increment_down(const double &deltaS) {
    down_counts++;
    down_dS += deltaS;
  }

  void print(FILE *fh) {

    fprintf(fh, "==== Statistics for %s ====\n", description.c_str());
    double tot = down_counts + up_counts;
    double prob = (up_yes + down_counts) / tot;
    fprintf(fh, "Probability of up acceptance = %f\n", prob);

    double mean_dS = up_dS / up_yes;
    prob = up_yes / tot;
    fprintf(fh,
            "dS >0: accepts = %ld; probability of up accept = %f; mean dS=%f\n",
            up_yes, prob, mean_dS);

    prob = up_no / tot;
    fprintf(fh,
            "dS =0: rejects = %ld; probability of up reject = %f; mean dS=%f\n",
            up_no, prob, 0.);

    prob = up_counts / tot;
    mean_dS = up_dS / up_counts;
    fprintf(fh,
            "dS>=0: counts  = %ld; probability of non-neg.  = %f; mean dS=%f\n",
            up_counts, prob, mean_dS);

    prob = down_counts / tot;
    mean_dS = down_dS / down_counts;
    fprintf(fh,
            "dS <0: counts  = %ld; probability of down steps= %f; mean dS=%f\n",
            down_counts, prob, mean_dS);

    fprintf(fh, "============================================\n");
  }
};

/////////////////////////////////////////////////////////////////////////
// Ideal Position Verlet Stepper with optional accept reject step
class IdealPV2 : public Stepper {
public:
  // If accept_reject is false then the step is done without
  // the metropolis accept reject
  IdealPV2(ModelA &in, const bool &accept_reject_in = false);
  void finalize() override;
  bool step(const double &dt) override;

  ~IdealPV2() { ; }

  PetscScalar computeEnergy(PetscScalar dt);
  void getEnergy(PetscScalar &oldE, PetscScalar &newE) const {
    oldE = oldEnergy;
    newE = newEnergy;
  }

private:
  ModelA *model;

  DM da;
  const ModelAData &data;

  bool step_no_reject(const double &dt);
  void rotatePhi(G_node ***phinew, double dt);

  bool accept_reject;
  o4_stepper_monitor monitor;
  PetscScalar oldEnergy;
  PetscScalar newEnergy;

  Vec previoussolution;

  // Minor optimization
  PetscInt xstart, ystart, zstart, xdimension, ydimension, zdimension;
};

/////////////////////////////////////////////////////////////////////////
// The heat bath step for the phi fields
class EulerLangevinHB : public Stepper {
public:
  EulerLangevinHB(ModelA &in);
  bool step(const double &dt);
  void finalize();
  ~EulerLangevinHB() { ; }

private:
  ModelA *model;
  o4_stepper_monitor monitor;

  Vec phi_local;
};

/////////////////////////////////////////////////////////////////////////
// The heat bath step for the charge fields
class ModelGChargeHB : public Stepper {
public:
  ModelGChargeHB(ModelA &in);
  bool step(const double &dt);
  void finalize();
  ~ModelGChargeHB() { ; }

private:
  ModelA *model;
  o4_stepper_monitor qmonitor;

  Vec phi_local;
  Vec dn_local;
};

// These are helper structures for ModelGChargeHB
struct g_face_case {
  int eoA; // zero or one / even or odd of site A
  int iB;  // either zero or one (indexing of the B cell relative to A)
  int jB;  // either zero or one  (indexing of the B cell relative to A)
  int kB;  // either zero or one  (indexing of the B cell relative to A)
  void print() {
    std::cout << "(EOA,iB,jB,kB): " << eoA << " (" << iB << "," << jB << ","
              << kB << ")" << std::endl;
  }
};

extern g_face_case g_face_cases[3][2];

PetscScalar modelg_update_charge_pair(const double &chi, const double &rms,
                                      const PetscScalar &nA,
                                      const PetscScalar &nB,
                                      o4_stepper_monitor &monitor);

/////////////////////////////////////////////////////////////////////////
class PV2HBSplit : public Stepper {
public:
  // The ideal step is A, heat bath step is B, the diffusion step is C.
  //
  // The format for  scounts = {inner, outer}.
  // With scounts  {1, 1} the stepping is (AB)C
  // With scounts  {3, 1} the stepping is (ABBB)C
  // With scounts  {1, 2} the stepping is (AB)(AB) C
  // With scounts  {2, 3} the stepping is (ABB)(ABB) (ABB) C
  //
  // One can switch off the diffusion step by setting nodiffuse=true
  //
  // One can only do diffusion by setting onlydiffuse=true
  PV2HBSplit(ModelA &in, const std::array<unsigned int, 2> &scounts = {1, 1},
             const bool &nodiffuse = false, const bool &onlydiffuse = false);
  bool step(const double &dt) override;
  void finalize() override {
    pv2.finalize();
    hbPhi.finalize();
    hbN.finalize();
  }
  const IdealPV2 &getIdealPV2() const { return pv2; }

  ~PV2HBSplit() { ; }

private:
  ModelA *model;

  IdealPV2 pv2;
  EulerLangevinHB hbPhi;
  ModelGChargeHB hbN;

  std::array<unsigned int, 2> stepcounts;
  bool nodiffusion;
  bool onlydiffusion;
};

#endif
