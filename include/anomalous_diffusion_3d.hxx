#pragma once
#ifndef ANOMALOUS_DIFFUSION_3D_H
#define ANOMALOUS_DIFFUSION_3D_H

#include "component.hxx"
#include "../include/div_ops.hxx"
#include <bout/fv_ops_impl.hxx>


/// Add anomalous diffusion of density, momentum and energy
///
/// # Mesh inputs
///
/// D_<name>, chi_<name>, nu_<name>
/// e.g `D_e`, `chi_e`, `nu_e`
///
/// in units of m^2/s
///
struct AnomalousDiffusion3D : public NamedComponent<AnomalousDiffusion3D> {
  /// # Inputs
  ///
  /// - <name>
  ///   - anomalous_D    This overrides D_<name> mesh input
  ///   - anomalous_chi  This overrides chi_<name>
  ///   - anomalous_nu   Overrides nu_<name>
  ///   - anomalous_sheath_flux  Allow anomalous flux into sheath?
  //                             Default false.
  AnomalousDiffusion3D(std::string name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "anomalous_diffusion_3d";

private:
  std::string name; ///< Species name

  bool diagnose;                           ///< Outputting diagnostics?
  bool include_D, include_D_par, include_chi, include_nu; ///< Which terms should be included?
  Field3D anomalous_D;                     ///< Anomalous density diffusion coefficient
  Field3D anomalous_chi;                   ///< Anomalous thermal diffusion coefficient
  Field3D anomalous_nu;                    ///< Anomalous momentum diffusion coefficient
  Field3D anomalous_D_par;                 ///< Anomalous parallel density diffusion coefficient
  std::shared_ptr<FCI::dagp_fv> dagp_op;
  
  bool anomalous_sheath_flux; ///< Allow anomalous diffusion into sheath?

  /// Inputs
  /// - species
  ///   - <name>
  ///     - AA
  ///     - density
  ///     - temperature  (optional)
  ///     - velocity     (optional)
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - <name>
  ///     - density_source
  ///     - momentum_source
  ///     - energy_source
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<AnomalousDiffusion3D> registercomponentanomalousdiffusion3d;
}

#endif // ANOMALOUS_DIFFUSION_3D_H
