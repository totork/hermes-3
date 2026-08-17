#pragma once
#ifndef LOWSOURCES_H
#define LOWSOURCES_H

#include "component.hxx"
#include "../include/div_ops.hxx"


struct LowSources : public NamedComponent<LowSources> {

  LowSources(std::string name, Options& alloptions, Solver*);

  static constexpr auto type = "low_sources";

private:
  std::string name; ///< Species name

  BoutReal low_density, low_density_width;
  BoutReal low_temperature, low_temperature_width;
  BoutReal low_timescale;
  
  void transform_impl(GuardedOptions& state) override;

  Field3D sourceterm(const Field3D ar, const BoutReal val, const BoutReal scale, const BoutReal width);
};

namespace {
RegisterComponent<LowSources> registercomponentlowsources;
}

#endif // LOWSOURCES_H
