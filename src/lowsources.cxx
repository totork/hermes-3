#include "../include/lowsources.hxx"

#include "../include/div_ops.hxx"

#include <bout/difops.hxx>
#include <bout/fv_ops.hxx>
#include <bout/fv_ops_impl.hxx>

#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/output_bout_types.hxx>

using bout::globals::mesh;

Field3D LowSources::sourceterm(const Field3D ar, const BoutReal val, const BoutReal scale,
                               const BoutReal width) {
  Field3D result = 0.0;

  const BoutReal eps = 1e-12;

  BOUT_FOR(i, ar.getRegion("RGN_NOY")) {
    // Compute ratio: 1 - (ar / val)
    BoutReal ratio = 1.0 - (ar[i] / val);

    // Smooth max(ratio, 0) using sqrt
    // This is: 0.5 * (x + sqrt(x² + width²)) for x = ratio
    BoutReal sqrt_term = std::sqrt(ratio * ratio + width * width);
    BoutReal smooth_max = 0.5 * (ratio + sqrt_term);

    result[i] = smooth_max / scale;
  }

  return result;
}

LowSources::LowSources(std::string name, Options& alloptions, Solver*)
    : NamedComponent(name, {readOnly("species:{name}:density", Regions::Interior),
                            readIfSet("species:{name}:{optional}", Regions::Interior),
                            readWrite("species:{name}:{output}")}) {
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];

  Options& options = alloptions[name];

  low_density = options["low_density"]
                    .doc("Density at which the low sources should activate?")
                    .withDefault<BoutReal>(-1.0)
                / Nnorm;

  low_density_width = options["low_density_width"]
                          .doc("Width of the transition for the low sources in density?")
                          .withDefault<BoutReal>(0.001);

  low_temperature = options["low_temperature"]
                        .doc("Temperature at which the low sources should activate?")
                        .withDefault<BoutReal>(-1.0)
                    / Tnorm;

  low_temperature_width =
      options["low_temperature_width"]
          .doc("Width of the transition for the low sources in temperature?")
          .withDefault<BoutReal>(0.001);

  low_timescale = options["low_timescale"]
                      .doc("Strength of the low source?")
                      .withDefault<BoutReal>(1e-6)
                  * Omega_ci;

  substitutePermissions("name", {name});
  substitutePermissions("optional", {"temperature", "velocity"});
  substitutePermissions("name", {name});
  substitutePermissions("optional", {"temperature", "velocity"});

  std::vector<std::string> output_vars;
  if (low_density > 0.0) {
    output_vars.push_back("density_source");
  }

  if (low_temperature > 0.0 || low_density > 0.0) {
    output_vars.push_back("energy_source");
  }

  substitutePermissions("output", output_vars);
}

void LowSources::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][objectName()];

  const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);

  const Field3D T = species.isSet("temperature")
                        ? GET_NOBOUNDARY(Field3D, species["temperature"])
                        : 0.0;
  if (low_density > 0.0) {
    Field3D N_src =
        LowSources::sourceterm(N, low_density, low_timescale, low_density_width);
    add(species["density_source"], N_src);
    if (species.isSet("pressure")) {
      add(species["energy_source"], (3.0 / 2.0) * T * N_src);
    }
  }

  if (low_temperature > 0.0) {
    Field3D E_src = (3.0 / 2.0) * N
                    * LowSources::sourceterm(T, low_temperature, low_timescale,
                                             low_temperature_width);
    add(species["energy_source"], E_src);
  }
}
