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
// Implements a deterministic diffusive step in Model G
class ModelGDiffusionStep : public Stepper {
public:
  ModelGDiffusionStep(ModelA &in, const bool &implicit_step = true);
  bool step(const double &dt);
  void finalize();
  ~ModelGDiffusionStep() { ; }
  static PetscErrorCode Form3PointLaplacian(DM da, Mat J, const double &hx,
                                            const double &hy, const double &hz,
                                            const double &GammaOverD);

private:
  ModelA *model;
  bool use_implicit_step; // If true, use implicit step, otherwise explicit

  Vec rhs;
  Vec dn;
  Mat J;
  Mat A;

  KSP ksp;
};

/////////////////////////////////////////////////////////////////////////
// Implements a deterministic explict diffusive step in Model G
class ModelGExplicitDiffusionStep : public Stepper {
public:
  ModelGExplicitDiffusionStep(ModelA &in) : model(&in) { ; }
  bool step(const double &dt) override;
  void finalize() override { ; }
  ~ModelGExplicitDiffusionStep() { ; }

private:
  ModelA *model;
};
/////////////////////////////////////////////////////////////////////////

// Solve the deterministic superfluid equations of motion using operator
// splitting
class SuperSplitStep : public Stepper {
public:
  SuperSplitStep(ModelA &in, const std::string &step_sequence = "AB",
                 const bool &use_implicit = false)
      : model(&in), pv2(in), steps(step_sequence), A_count(0), B_count(0) {
    if (use_implicit) {
      diffusion = std::make_unique<ModelGDiffusionStep>(in, true);
    } else {
      diffusion = std::make_unique<ModelGExplicitDiffusionStep>(in);
    }
    for (char s : steps) {
      if (s == 'A') {
        A_count++;
      } else if (s == 'B') {
        B_count++;
      }
    }
  }
  bool step(const double &dt) override {
    PetscLogEvent ideal_log, diffusive_log;
    PetscLogEventRegister("ExplicitDiffusion", 0, &diffusive_log);
    PetscLogEventRegister("IdealStep", 0, &ideal_log);

    bool ok = true;
    for (char s : steps) {
      if (s == 'A') {
        PetscLogEventBegin(ideal_log, 0, 0, 0, 0);
        ok = ok && pv2.step(dt / A_count);
        PetscLogEventEnd(ideal_log, 0, 0, 0, 0);
      } else if (s == 'B') {
        PetscLogEventBegin(diffusive_log, 0, 0, 0, 0);
        ok = ok && diffusion->step(dt / B_count);
        PetscLogEventEnd(diffusive_log, 0, 0, 0, 0);
      }
    }
    return ok;
  }
  void finalize() override {
    pv2.finalize();
    diffusion->finalize();
  }

private:
  ModelA *model;
  IdealPV2 pv2;
  std::unique_ptr<Stepper> diffusion;
  std::string steps;  // The steps to take, e.g. "ABBB"
  double A_count = 0; // Number of A steps in the string steps
  double B_count = 0; // Number of B steps in the string steps
};

/////////////////////////////////////////////////////////////////////////
class PV2HBSplit : public Stepper {
public:
  // The ideal step is A, heat bath step is B, the diffusion step is C.
  //
  // The order of the steps is specified with a string, e.g. "ABC"
  // takes a step of A, then B, then C  each for size dt
  //
  // The string which was primarily used in production model g code is
  // "ABBABBABBC" . In this case, the timestep for B = dt/6, the timestep for A
  // = dt/3, and the timestep for C = dt.
  //
  // The default is to include all steps, i.e. ideal, heat bath and diffusion.
  // But the individual steps can be turned off by setting the flags, ideal,
  // heatbath, and diffusion to false. This can be done at construction
  // PV2HBSplit(ModelA &in, "ABBABBABBC", true, true, true)  or later with
  // setmode().
  PV2HBSplit(ModelA &in, const std::string &inputsteps,
             const bool &ideal = true, const bool &heatbath = true,
             const bool &diffusion = true);
  bool step(const double &dt) override;
  void finalize() override {
    pv2.finalize();
    hbPhi.finalize();
    hbN.finalize();
  }
  const IdealPV2 &getIdealPV2() const { return pv2; }

  // Set the mode of the stepper. If ideal is false, then the ideal step is not
  // performed, etc. The default, which is all true, is set with setmode().
  void setmode(const bool &ideal = true, const bool &heatbath = true,
               const bool &diffusion = true) {
    include_ideal = ideal;
    include_heatbath = heatbath;
    include_diffusion = diffusion;
  }

  ~PV2HBSplit() { ; }

private:
  ModelA *model;

  IdealPV2 pv2;
  EulerLangevinHB hbPhi;
  ModelGChargeHB hbN;

  std::string steps; // The steps to take, e.g. "ABC"

  bool include_ideal;     // If true, do the ideal step
  bool include_heatbath;  // If true, do the heat bath step
  bool include_diffusion; // If true, do the diffusion step
  double A_count;         // Number of A steps in the string steps
  double B_count;         // Number of B steps in the string steps
  double C_count;         // Number of C steps in the string steps
};

#endif
