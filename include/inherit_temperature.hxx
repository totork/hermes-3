#pragma once
#ifndef INHERIT_TEMPERATURE
#define INHERIT_TEMPERATURE

#include "component.hxx"

struct Inherittemperature : public Component {
  /// Inputs
  /// ------
  ///
  /// @param name     Short name for species e.g. "e"
  /// @param alloptions  Component configuration options
  ///   - <name>
  ///     - charge   Required to have a particle charge
  ///     - AA       Atomic mass
  ///
  Inherittemperature(std::string name, Options &alloptions, Solver *UNUSED(solver));

  /// 
  /// Sets in state
  /// - species
  ///   - <name>
  ///     - density
  ///     - charge
  ///     - AA
  void transform(Options &state) override;

  /// Get the final density for output
  /// including any boundary conditions applied
  void finally(const Options &state) override;

  void outputVars(Options &state) override;
private:
  std::string name; ///< Name of this species
  std::string donor_species;
  Field3D temp;
};

namespace {
RegisterComponent<Inherittemperature>
    registercomponentquasineutral("inherit_temperature");
}

#endif // INHERIT_TEMPERATURE
