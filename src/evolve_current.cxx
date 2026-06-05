
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/evolve_current.hxx"
#include "../include/div_ops.hxx"
#include "../include/hermes_build_config.hxx"


using bout::globals::mesh;

EvolveCurrent::EvolveCurrent(std::string name, Options &alloptions, Solver *solver) {
  AUTO_TRACE();
  
  // Evolve the momentum in time
  solver->add(VePsi, std::string("VePsi"));

  auto& options = alloptions[name];

  exb_advection = options["exb_advection"]
                   .doc("Include ExB advection?")
                   .withDefault<bool>(true);


  const Options& units = alloptions["units"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Lnorm = units["meters"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];

  
  auto Cs0 = Lnorm * Omega_ci;
  BoutReal lambda_ei = 24. - log(sqrt(Nnorm / 1e6) / Tnorm);
  BoutReal tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * lambda_ei * pow(Tnorm, -3. / 2));
  tau_e1 = (Cs0 / Lnorm ) * tau_e0;

  

  diagnose = options["diagnose"]
    .doc("Output additional diagnostics?")
    .withDefault<bool>(false);

  main_ionspecies = options["main_ionspecies"]
    .doc("What is the main ion species the current should model the displacement from?")
    .withDefault<std::string>("h+");

  viscosity = options["viscosity"]
    .doc("Use parallel viscosity for velocity?")
    .withDefault<bool>(false);

  

  if (mesh->isFci()) {
    const auto coord = mesh->getCoordinates();
    // Note: This is 1 for a Clebsch coordinate system
    //       Remove parallel slices before operations
    bracket_factor = sqrt(coord->g_22.withoutParallelSlices())
      / (coord->J.withoutParallelSlices() * coord->Bxy);
  } else {
    // Clebsch coordinate system
    bracket_factor = 1.0;
  }

}

void EvolveCurrent::transform(Options &state) {
  AUTO_TRACE();
  auto& species = state["species"]["e"];
  auto& main_ions = state["species"][main_ionspecies];
  BoutReal AA = get<BoutReal>(species["AA"]);


  VePsi.applyBoundary();
  mesh->communicate(VePsi);
  VePsi.applyParallelBoundary();
  
  auto N = getNoBoundary<Field3D>(species["density"]);
  Field3D Vi = getNoBoundary<Field3D>(main_ions["velocity"]);
  Ve = VePsi + Vi;

  set(species["velocity"], Ve);
  set(species["momentum"], Ve*N*AA);

  
}

void EvolveCurrent::finally(const Options &state) {
  AUTO_TRACE();

  auto& species = state["species"]["e"];
  auto& main_ions = state["species"][main_ionspecies];
  
  BoutReal AA = get<BoutReal>(species["AA"]);
  BoutReal AA_ions = get<BoutReal>(main_ions["AA"]);

  BoutReal mi_me = (AA_ions / AA);
  
  Field3D T = get<Field3D>(species["temperature"]);
  Field3D N = get<Field3D>(species["density"]);
  Field3D P = get<Field3D>(species["pressure"]);

  Field3D Vi = get<Field3D>(main_ions["velocity"]);
  Ve = get<Field3D>(species["velocity"]);
  
  ddt(VePsi) = 0.0;

  if (state["fields"].isSet("phi")) {
    const Field3D phi = get<Field3D>(state["fields"]["phi"]);
    ddt(VePsi) += mi_me * Grad_par(phi);
  }

  ddt(VePsi) -= 0.71 * mi_me * Grad_par(T);

  ddt(VePsi) -= mi_me * Grad_par(P) / N;

  // Use Spitzer resistivity here
  Field3D T32 = T * sqrt(T);
  Field3D tau_e = tau_e1 * T32 / N;
  Field3D nu = 1.0 / (1.96 * tau_e * mi_me);

  ddt(VePsi) += mi_me * nu * (Vi - Ve);


  if (viscosity) {

    Field3D eta_epar = 0.973 * mi_me * tau_e * T;
    eta_epar.applyBoundary("neumann");
    mesh->communicate(eta_epar);
    eta_epar.applyParallelBoundary("parallel_neumann_o1");

    Field3D dummy;
    ddt(VePsi) += Div_par_K_Grad_par_mod(eta_epar,Ve, dummy, true);
    
  }
  
  
  
}

void EvolveCurrent::outputVars(Options &state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  set_with_attrs(state[std::string("Ve")], Ve,
                   {{"time_dimension", "t"},
                    {"units", "m / s"},
                    {"conversion", Cs0},
                    {"long_name", " parallel electron velocity"},
                    {"standard_name", "electron velocity"},
                    {"source", "evolve_current"}});
  
}
