
#include <numeric> // for accumulate

#include "../include/inherit_temperature.hxx"

Inherittemperature::Inherittemperature(std::string name, Options &alloptions,
                           Solver *UNUSED(solver))
    : name(name) {
  Options &options = alloptions[name];

  donor_species = options["donor_T"].withDefault("h+");
}

void Inherittemperature::transform(Options &state) {
  AUTO_TRACE();
  // Iterate through all subsections
  Options &allspecies = state["species"];
  Options &donor = allspecies[donor_species];
  Options &species = allspecies[name];

  
  Field3D donor_T = getNoBoundary<Field3D>(donor["temperature"]);
  
  set(species["temperature"], donor_T);

  Field3D den = getNoBoundary<Field3D>(species["density"]);
  
  set(species["pressure"], den + donor_T);
  
}

void Inherittemperature::finally(const Options &state) {
  temp = get<Field3D>(state["species"][name]["temperature"]);
}

void Inherittemperature::outputVars(Options &state) {
  AUTO_TRACE();
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);

  set_with_attrs(state[std::string("T") + name], temp,
                 {{"time_dimension", "t"},
                  {"units", "eV"},
                  {"conversion", Tnorm},
                  {"source", "inherit_temperature"}});
  
  
}
