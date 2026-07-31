#include "../include/anomalous_diffusion_3d.hxx"

#include "../include/div_ops.hxx"

#include <bout/fv_ops_impl.hxx>
#include <bout/fv_ops.hxx>
#include <bout/difops.hxx>

#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/output_bout_types.hxx>

using bout::globals::mesh;

AnomalousDiffusion3D::AnomalousDiffusion3D(std::string name, Options& alloptions, Solver*)
    : NamedComponent(name, {readOnly("species:{name}:density", Regions::Interior),
                            readIfSet("species:{name}:{optional}", Regions::Interior),
                            readWrite("species:{name}:{output}")}) {
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  const BoutReal diffusion_norm = rho_s0 * rho_s0 * Omega_ci; // m^2/s

  Options& options = alloptions[name];

  // Set in the mesh or options (or both)
  anomalous_D = 0.0;
  include_D = (mesh->get(anomalous_D, std::string("D_") + name) == 0)
              || options.isSet("anomalous_D");
  // Option overrides mesh value
  anomalous_D = options["anomalous_D"]
                    .doc("Anomalous particle diffusion coefficient [m^2/s]")
                    .withDefault(anomalous_D)
                / diffusion_norm;

  anomalous_D_par = 0.0;
  include_D_par = (mesh->get(anomalous_D_par, std::string("D_") + name) == 0)
              || options.isSet("anomalous_D_par");
  // Option overrides mesh value 
  anomalous_D_par = options["anomalous_D_par"]
                    .doc("Anomalous parallel particle diffusion coefficient [m^2/s]")
                    .withDefault(anomalous_D_par)
                / diffusion_norm;

  if (mesh->isFci()) {
    mesh->communicate(anomalous_D_par);
    anomalous_D_par.applyParallelBoundary("parallel_neumann_o1");
  }

  anomalous_chi = 0.0;
  include_chi = (mesh->get(anomalous_chi, std::string("chi_") + name) == 0)
                || options.isSet("anomalous_chi");
  anomalous_chi = options["anomalous_chi"]
                      .doc("Anomalous thermal diffusion coefficient [m^2/s]")
                      .withDefault(anomalous_chi)
                  / diffusion_norm;

  anomalous_nu = 0.0;
  include_nu = (mesh->get(anomalous_nu, std::string("nu_") + name) == 0)
               || options.isSet("anomalous_nu");
  anomalous_nu = options["anomalous_nu"]
                     .doc("Anomalous momentum diffusion coefficient [m^2/s]")
                     .withDefault(anomalous_nu)
                 / diffusion_norm;

  anomalous_sheath_flux = options["anomalous_sheath_flux"]
                              .doc("Allow anomalous diffusion into sheath?")
                              .withDefault<bool>(false);

  diagnose = alloptions[name]["diagnose"]
                 .doc("Output additional diagnostics?")
                 .withDefault<bool>(false);

#if BOUT_USE_METRIC_3D
  dagp_op = FCI::getDagp_fv(mesh, rho_s0);
#endif 
  
  substitutePermissions("name", {name});
  substitutePermissions("optional", {"temperature", "velocity"});
  std::vector<std::string> output_vars;
  if (include_D || include_D_par) {
    output_vars.push_back("density_source");
    output_vars.push_back("particle_flow_xlow");
    output_vars.push_back("particle_flow_ylow");
  }
  if (include_D or include_chi) {
    output_vars.push_back("energy_source");
    output_vars.push_back("energy_flow_xlow");
    output_vars.push_back("energy_flow_ylow");
  }
  if (include_D or include_nu) {
    setPermissions(readOnly(fmt::format("species:{}:AA", name)));
    output_vars.push_back("momentum_source");
    output_vars.push_back("momentum_flow_xlow");
    output_vars.push_back("momentum_flow_ylow");
  }
  substitutePermissions("output", output_vars);
}

void AnomalousDiffusion3D::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][objectName()];

  // Diffusion operates on 2D (axisymmetric) profiles
  // Note: Includes diffusion in Y, so set boundary fluxes
  // to zero by imposing neumann boundary conditions.
  const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
  Field2D N2D = DC(N);

  const Field3D T = species.isSet("temperature")
                        ? GET_NOBOUNDARY(Field3D, species["temperature"])
                        : 0.0;


  const Field3D V =
      species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"]) : 0.0;


  Field3D flow_xlow, flow_ylow; // Flows through cell faces

  if (include_D) {
    // Particle diffusion. Gradients of density drive flows of particles,
    // momentum and energy. The implementation here is equivalent to an
    // advection velocity
    //
    //  v_D = - D Grad_perp(N) / N

    add(species["density_source"],
        (*dagp_op)(anomalous_D, N, flow_xlow, flow_ylow, false));
    add(species["particle_flow_xlow"], flow_xlow);
    add(species["particle_flow_ylow"], flow_ylow);

    // Note: Upwind operators used, or unphysical increases
    // in temperature and flow can be produced
    auto AA = get<BoutReal>(species["AA"]);
    add(species["momentum_source"],
        (*dagp_op)(AA * V * anomalous_D, N, flow_xlow,
		flow_ylow, false));
    add(species["momentum_flow_xlow"], flow_xlow);
    add(species["momentum_flow_ylow"], flow_ylow);

    add(species["energy_source"],
        (*dagp_op)((3. / 2) * T * anomalous_D, N,
		flow_xlow, flow_ylow, false));
    add(species["energy_flow_xlow"], flow_xlow);
    add(species["energy_flow_ylow"], flow_ylow);
  }

  if (include_D_par) {
    Field3D dummy;
    add(species["density_source"], Div_par_K_Grad_par_mod(anomalous_D_par, N, dummy, false));
  }

  if (include_chi) {
    // Gradients in temperature that drive energy flows
    add(species["energy_source"],
        (*dagp_op)(anomalous_chi * N, T, flow_xlow,
		   flow_ylow, false));
    add(species["energy_flow_xlow"], flow_xlow);
    add(species["energy_flow_ylow"], flow_ylow);
  }

  if (include_nu) {
    // Gradients in flow speed that drive momentum flows
    auto AA = get<BoutReal>(species["AA"]);
    add(species["momentum_source"],
        (*dagp_op)(anomalous_nu * AA * N, V, flow_xlow,
		flow_ylow, false));
    add(species["momentum_flow_xlow"], flow_xlow);
    add(species["momentum_flow_ylow"], flow_ylow);
  }
}

void AnomalousDiffusion3D::outputVars(Options& state) {
  // Normalisations
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  const std::string& name = objectName();

  if (diagnose) {

    // Save particle, momentum and energy channels

    set_with_attrs(state[{std::string("anomalous_D_") + name}], anomalous_D,
                   {{"time_dimension", "t"},
                    {"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "anomalous density diffusion"},
                    {"long_name", std::string("Anomalous density diffusion of ") + name},
                    {"source", "anomalous_diffusion"}});

    set_with_attrs(state[{std::string("anomalous_Chi_") + name}], anomalous_chi,
                   {{"time_dimension", "t"},
                    {"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "anomalous thermal diffusion"},
                    {"long_name", std::string("Anomalous thermal diffusion of ") + name},
                    {"source", "anomalous_diffusion"}});

    set_with_attrs(state[{std::string("anomalous_nu_") + name}], anomalous_nu,
                   {{"time_dimension", "t"},
                    {"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "anomalous momentum diffusion"},
                    {"long_name", std::string("Anomalous momentum diffusion of ") + name},
                    {"source", "anomalous_diffusion"}});
  }
}
