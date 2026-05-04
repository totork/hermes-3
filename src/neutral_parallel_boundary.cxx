#include "../include/neutral_parallel_boundary.hxx"


#include <bout/output_bout_types.hxx>

#include "bout/constants.hxx"
#include "bout/mesh.hxx"

#include "bout/parallel_boundary_region.hxx"
#include "bout/boundary_iterator.hxx"

using bout::globals::mesh;


NeutralParallelBoundary::NeutralParallelBoundary(std::string name, Options &alloptions, Solver *) {
  AUTO_TRACE();
  
  Options &options = alloptions[name];
  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];

  boundary_value_yup = options["boundary_value_yup"]
           .doc("Boundary value for the neutral density")
           .withDefault(-1.0) / Nnorm;

  boundary_value_ydown = options["boundary_value_ydown"]
           .doc("Boundary value for the neutral density")
           .withDefault(-1.0) / Nnorm;

  yboundary.init(options);
  
}



void NeutralParallelBoundary::transform(Options &state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];


  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if ((species.isSet("charge") ))
      continue; // Skip when species has charge

    Field3D N = toFieldAligned(floor(GET_NOBOUNDARY(Field3D, species["density"]), 0.0));

    iter_regions([&](auto& region) {
      for (auto& pnt : region) {
	const auto& i = pnt.ind();
	if (pnt.dir > 0.0) {
	  pnt.dirichlet_o2(N, boundary_value_yup);
	} else {
	  pnt.dirichlet_o2(N, boundary_value_ydown);
	}

      }
    });
  
    setBoundary(species["density"], fromFieldAligned(N));
    
  }
  
}
