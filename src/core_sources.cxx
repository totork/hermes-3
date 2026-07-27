/// Braginskii electron viscosity

#include <bout/fv_ops.hxx>
#include <bout/mesh.hxx>
#include <bout/difops.hxx>
#include <bout/constants.hxx>
#include "../include/div_ops.hxx"

#include "../include/core_sources.hxx"


using bout::globals::mesh;


CoreSources::CoreSources(std::string name, Options& alloptions, Solver*) {
  auto& options = alloptions[name];

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Lnorm = units["meters"];
  const BoutReal Anorm = (Lnorm * Lnorm);

  BoutReal eV = 6.241509e+18;

  const BoutReal Particlenorm = Nnorm * Omega_ci;
  const BoutReal Powernorm = Nnorm * Tnorm * Omega_ci / eV * (Lnorm * Lnorm * Lnorm);
  
  // diagnose = options["diagnose"].doc("Output diagnostics?").withDefault<bool>(false);


  input_power = options["input_power"].doc("Input power for each species? [MW]").withDefault<BoutReal>(-1.0) / Powernorm;
  input_particleflux = options["input_particleflux"].doc("Input particleflux for each species?").withDefault<BoutReal>(-1.0) / Particlenorm;  
  core_area = options["core_area"].doc("Core area that is used to divide the total input power to get the power density?").withDefault<BoutReal>(-1.0) / Anorm;

  
  
  densitytarget = options["densitytarget"].doc("Target density the sources should be cut off after?").withDefault<BoutReal>(-1.0) / Nnorm;
  temperaturetarget = options["temperaturetarget"].doc("Target density the sources should be cut off after?").withDefault<BoutReal>(-1.0)/ Tnorm;
  target_timescale = options["target_timescale"].withDefault<BoutReal>(0.02);

  BoutReal powerfluxdensity = input_power / core_area;
  BoutReal particlefluxdensity = input_particleflux / core_area;
  
  source_density = 0.0;
  source_pressure = 0.0;


  Coordinates* coord = mesh->getCoordinates();

  Field3D cellcross_x = coord->cell_volume() / (coord->dx * sqrt(coord->g_11));
  Field3D cellcross_z = coord->cell_volume() / (coord->dz * sqrt(coord->g_33));
  
  if (mesh->firstX()) {
    for (int i = mesh->xstart - 2; i >= 0; --i) {
      for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        for (int k = 0; k < mesh->LocalNz; ++k) {
	  BoutReal cellarea_xdown = 0.5 * (cellcross_x(mesh->xstart, j, k) + cellcross_x(mesh->xstart, j, k));
	  
	  source_density(mesh->xstart, j, k) = particlefluxdensity * cellarea_xdown / coord->cell_volume()(mesh->xstart, j, k);
	  source_pressure(mesh->xstart, j, k) = powerfluxdensity * cellarea_xdown / coord->cell_volume()(mesh->xstart, j, k);
        }
      }
    }
  }
  
  
  
}



void CoreSources::transform(Options& state) {

  Options& allspecies = state["species"];
  
  for (auto& kv : allspecies.getChildren()) {
    
    const auto& species_name = kv.first;

    Options& species = allspecies[species_name];


    if (std::fabs(get<BoutReal>(species["charge"])) < 1e-3) {
      // No charge                                                                                                                                                                                                
      continue;
    }

    if (std::fabs(get<BoutReal>(species["charge"])) > 1.01) {
      // Not hydrogen / hydrogen isotope                                                                                                                                                                          
      continue;
    }

    Field3D species_density_source = species.isSet("density_source")
                                 ? getNonFinal<Field3D>(species["density_source"])
                                 : 0.0;

    Field3D species_energy_source = species.isSet("energy_source")
                                 ? getNonFinal<Field3D>(species["energy_source"])
                                 : 0.0;

    const Field3D N = get<Field3D>(species["density"]);
    const Field3D P = get<Field3D>(species["pressure"]);
    const Field3D T = P / N;
    
    if (densitytarget > 0.0) {
      species_density_source = species_density_source + adaptive_sourceterm(N, source_density, densitytarget, target_timescale);
    } else {
      species_density_source = species_density_source + source_density;
    }


    if (temperaturetarget > 0.0) {
      species_energy_source = species_energy_source + adaptive_sourceterm(T, source_pressure, temperaturetarget, target_timescale);
    } else {
      species_energy_source = species_energy_source + source_pressure;
    }
       
    
    set<Field3D>(species["density_source"], species_density_source);
    set<Field3D>(species["energy_source"], species_energy_source);
    
  }
  
}

void CoreSources::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  set_with_attrs(state["core_source_density"], source_density,
		 {{"source", "core_sources"}});
  
  set_with_attrs(state["core_source_pressure"], source_pressure,
		 {{"source", "core_sources"}});
  



  
}
