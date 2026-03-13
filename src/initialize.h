#ifndef INITIALIZE_H
#define INITIALIZE_H
#include "ModelA.h"

void initialize_gaussians(G_node *node, const double &x, const double &y,
                          const double &z, ModelA *model, void *ctx);
void initialize_wave_spins(G_node *node, const double &x, const double &y,
                           const double &z, ModelA *model, void *ctx);
#endif
