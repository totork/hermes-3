
#pragma once
#ifndef EVOLVE_CURRENT_H
#define EVOLVE_CURRENT_H

#include "component.hxx"

/// Evolve parallel momentum
struct EvolveCurrent : public Component {
  EvolveCurrent(std::string name, Options &options, Solver *solver);

  /// This sets in the state
  /// - species
  ///   - <name>
  ///     - momentum
  ///     - velocity  if density is defined
  void transform(Options &state) override;
  
  /// Calculate ddt(NV).
  ///
  /// Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///     - velocity
  ///     - pressure (optional)
  ///     - momentum_source (optional)
  ///     - sound_speed (optional, used for numerical dissipation)
  ///     - temperature (only needed if sound_speed not provided)
  /// - fields
  ///   - phi (optional)
  void finally(const Options &state) override;

  void outputVars(Options &state) override;
private:
  std::string name;     ///< Short name of species e.g "e"
  std::string Vname;
  Field3D VePsi;
  Field3D Ve;            ///< Species parallel velocity
  std::string main_ionspecies;
  bool diagnose;
  bool viscosity;
  bool exb_advection;
  Coordinates::FieldMetric bracket_factor; ///< For non-Clebsch coordinate systems (e.g. FCI)
  BoutReal tau_e1;
};

namespace {
RegisterComponent<EvolveCurrent> registercomponentevolvemomentum("evolve_current");
}

#endif // EVOLVE_CURRENT_H
