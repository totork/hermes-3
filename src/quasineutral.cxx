
#include <bout/mesh.hxx>

#include "../include/quasineutral.hxx"

#include <numeric> // for accumulate
#include <string>
using bout::globals::mesh;

Quasineutral::Quasineutral(std::string name, Options& alloptions, Solver* UNUSED(solver))
    : NamedComponent(name,
                     {readWrite("species:{name}:{outputs}"),
                      // FIXME: These are only read if BOTH are set
                      readIfSet("species:{all_species}:charge"),
                      readIfSet("species:{all_species}:density", Regions::Interior)}),
      name(name) {
  Options& options = alloptions[name];

  // Need to have a charge and mass
  charge = options["charge"].doc("Particle charge. electrons = -1");
  AA = options["AA"].doc("Particle atomic mass. Proton = 1");

  ASSERT0(charge != 0.0);
  substitutePermissions("name", {name});
  substitutePermissions("outputs", {"AA", "charge", "density"});
}

void Quasineutral::transform_impl(GuardedOptions& state) {
  // Iterate through all subsections
  GuardedOptions allspecies = state["species"];
  std::map<std::string, GuardedOptions> children = allspecies.getChildren();

  // Add charge density of other species
  Field3D rho;

  if (mesh->isFci()) {
    rho = 0.0;
    rho.applyBoundary("neumann");
    mesh->communicate(rho);
    rho.applyParallelBoundary("parallel_neumann_o1");
    for (auto& kv : children) {
      GuardedOptions species = allspecies[kv.first]; // Note: Need non-const
      
      if (kv.first == name) {
	continue;
      }

      if (!species.isSet("charge")) {
	continue;
      }
      
      auto q = get<BoutReal>(species["charge"]);
      if (fabs(q) < 1e-5 ) {
	continue;
      }

      if (!species.isSet("density")) {
	continue;
      }
      Field3D new_rho = getNoBoundary<Field3D>(species["density"]).asField3DParallel() * q;
      rho = rho.asField3DParallel() + new_rho.asField3DParallel();
      
    }
    ASSERT2(rho.hasParallelSlices());
    
  } else {
    rho= std::accumulate(
			 // Iterate through species
			 begin(children), end(children),
			 // Start with no charge
			 Field3D(0.0),
			 [this](Field3D value,
				const std::map<std::string, GuardedOptions>::value_type& name_species) {
			   const GuardedOptions species = name_species.second;
			   // Add other species which have density and charge
			   if (name_species.first != name and species.isSet("charge")
			       and species.isSet("density")) {
			     // Note: Not assuming that the boundary has been set
			     return Field3D{value
			       + getNoBoundary<Field3D>(species["density"])
                               * get<BoutReal>(species["charge"])};
			   }
			   return value;
			 });
  }

  // Set quantites for this species
  GuardedOptions species = allspecies[name];

  // Calculate density required. Floor so that density is >= 0
  if (mesh->isFci()) {
    density = floor(rho.asField3DParallel() / (-charge), 0.0);
  } else {
    density = floor(Field3D{rho / (-charge)}, 0.0);
  }
  
  set(species["density"], density);

  set(species["charge"], charge);

  set(species["AA"], AA);
}

void Quasineutral::finally(const Options& state) {
  // Density may have had boundary conditions applied
  density = get<Field3D>(state["species"][name]["density"]);
}

void Quasineutral::outputVars(Options& state) {
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
