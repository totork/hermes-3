#pragma once
#ifndef CORE_SOURCES_H
#define CORE_SOURCES_H

#include "component.hxx"

/// Electron viscosity 
///
/// Adds Braginskii parallel electron viscosity, with SOLPS-style
/// viscosity flux limiter
///
/// Needs to be calculated after collisions, because collision
/// frequency is used to calculate parallel viscosity
///
/// References
///  - https://farside.ph.utexas.edu/teaching/plasma/lectures1/node35.html
///
struct CoreSources : public Component {
  /// Inputs
  /// - <name>
  ///   - diagnose: bool, default false
  ///     Output diagnostic SNVe_viscosity?
  ///   - eta_limit_alpha: float, default -1.0
  ///     Flux limiter coefficient. < 0 means no limiter
  CoreSources(std::string name, Options& alloptions, Solver*);

  /// Inputs
  /// - species
  ///   - e
  ///     - pressure  (skips if not present)
  ///     - velocity  (skips if not present)
  ///     - collision_frequency
  ///
  /// Sets in the state
  /// - species
  ///   - e
  ///     - momentum_source
  ///
  void transform(Options &state) override;

  void outputVars(Options &state) override;
private:
  BoutReal input_power, input_particleflux, core_area;
  Field3D inner_area;
  BoutReal densitytarget, temperaturetarget;
  Field3D source_pressure, source_density;
};

namespace {
RegisterComponent<CoreSources> registercomponentelectronviscosity("core_sources");
}

#endif
