#include "ModelA.h"
#include <cmath>
#include <tuple>

// Rotate around Z-axis by phi (azimuthal), then around Y-axis by theta (polar)
std::tuple<double, double, double> rotateVector(double dx, double dy, double dz,
                                                double theta, double phi) {
  // First rotation: around Z-axis by phi (azimuthal angle)
  double cos_phi = std::cos(phi);
  double sin_phi = std::sin(phi);

  double x1 = dx * cos_phi - dy * sin_phi;
  double y1 = dx * sin_phi + dy * cos_phi;
  double z1 = dz;

  // Second rotation: around Y-axis by theta (polar angle)
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);

  double x2 = x1 * cos_theta + z1 * sin_theta;
  double y2 = y1;
  double z2 = -x1 * sin_theta + z1 * cos_theta;

  return std::make_tuple(x2, y2, z2);
}

void initialize_gaussians(G_node *node, const double &x, const double &y,
                          const double &z, ModelA *model, void *ctx) {
  auto &inputs = *reinterpret_cast<nlohmann::json *>(ctx);
  double sigmax = inputs["sigmax"];
  double sigmay = inputs["sigmay"];
  double sigmaz = inputs["sigmaz"];
  double A = inputs["amplitude"];
  data_node *u = reinterpret_cast<data_node *>(node);

  for (int L = 0; L < ModelAData::Ndof; L++) {
    double dx = x - model->data.LX / 2.0;
    double dy = y - model->data.LY / 2.0;
    double dz = z - model->data.LZ / 2.0;
    double theta = inputs["theta"];
    double phi = inputs["phi"];
    auto [xp, yp, zp] = rotateVector(dx, dy, dz, theta, phi);
    u->x[L] = A * exp(-xp * xp / (2. * sigmax * sigmax) -
                      yp * yp / (2. * sigmay * sigmay) -
                      zp * zp / (2. * sigmaz * sigmaz));
  }
}
// Routine that initializes the fields according to a wave with wave number k =
// n * 2 Pi/L and normalization phi^2 = f
//
// The specific form of the wave depends on the test_case parameter, which is
// passed through the input file. The inputfile is the "context" for this
// routine.
void initialize_wave_spins(G_node *node, const double &x, const double &y,
                           const double &z, ModelA *model, void *ctx) {
  auto &inputs = *reinterpret_cast<nlohmann::json *>(ctx);
  PetscScalar *u = reinterpret_cast<PetscScalar *>(node);

  int test_case = inputs["test_case"];
  double f = sqrt(model->data.f2());

  int wave_number = inputs.value("wave_number", 2);
  PetscReal k = wave_number * 2 * M_PI / model->data.LX;

  // PetscScalar chi = model->data.acoefficients.chi;
  for (int L = 0; L < ModelAData::Ndof; L++) {
    u[L] = 0.0; // Initialize all components to zero
  }

  if (test_case == 1) {
    u[0] = f * cos(k * x);
    u[1] = f * sin(k * x);
    // u[4] = k * chi ; // n_01
  } else if (test_case == 2) {
    u[0] = f * cos(k * x) * cos(k * y);
    u[1] = f * sin(k * x) * cos(k * y);
    u[2] = f * sin(k * y);
    // u[4] = k * chi ; // n_01
    // u[5] = k * chi ; // n_02
    // u[9] = k * chi ; // n_12 or nv[2]
  } else if (test_case == 3) {
    double a = M_PI;
    u[0] = f * cos(a * sin(k * x));
    u[1] = f * sin(a * sin(k * x));
  } else if (test_case == 4) {
    // Sod problem
    if (x >= model->data.LX / 4. and x < 3 * model->data.LX / 4.) {
      u[0] = f;
    } else {
      u[1] = f;
    }
  } else {
    throw std::runtime_error("Unknown test_case in initialize_wave_spins");
  }
}
