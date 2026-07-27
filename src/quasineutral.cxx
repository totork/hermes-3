
#include <numeric> // for accumulate

#include "../include/quasineutral.hxx"

Quasineutral::Quasineutral(std::string name, Options &alloptions,
                           Solver *UNUSED(solver))
    : name(name) {
  Options &options = alloptions[name];

  // Need to have a charge and mass
  charge = options["charge"].doc("Particle charge. electrons = -1");
  AA = options["AA"].doc("Particle atomic mass. Proton = 1");

  ASSERT0(charge != 0.0);
}

void Quasineutral::transform(Options &state) {
  // Iterate through all subsections
  Options &allspecies = state["species"];

  Field3D rho = 0.0;
  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (kv.first == name) {
      continue;
    }

    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5 ) {
      continue;
    }

    if (!species.isSet("density")) {
      continue;
    }

    rho = rho + getNoBoundary<Field3D>(species["density"]) * q;
    
  }

    

  // Set quantites for this species
  Options &species = allspecies[name];

  // Calculate density required. Floor so that density is >= 0
  Field3D den = rho / (-charge);
  density = floor(den, 0.0);
  set(species["density"], density);

  set(species["charge"], charge);

  set(species["AA"], AA);
}

void Quasineutral::finally(const Options &state) {
  // Density may have had boundary conditions applied
  density = get<Field3D>(state["species"][name]["density"]);
}

void Quasineutral::outputVars(Options &state) {
  auto Nnorm = get<BoutReal>(state["Nnorm"]);

  // Save the density
  set_with_attrs(state[std::string("N") + name], density,
                 {{"time_dimension", "t"},
                  {"units", "m^-3"},
                  {"conversion", Nnorm},
                  {"long_name", name + " number density"},
                  {"standard_name", "density"},
                  {"species", name},
                  {"source", "quasineutral"}});
}
