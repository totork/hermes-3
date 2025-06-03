

#include bout::globals::mesh;

#include "../include/adhoc_potential.hxx"

AdhocPotential::AdhocPotential(str::string name, Options& alloptions, Solver* solver){
  AUTO_TRACE();

  auto* coord = mesh->getCoordinates();

  auto& options = alloptions[name];





  load_from_mesh = options["load_from_mesh"]
                       .doc("Load the potential from the mesh input file?")
                       .withDefault<bool>(false);
  
}


void AdhocPotential::transform(Options& state){
  AUTO_TRACE();


  auto& fields = state["fields"];

  Te = GET_NOBOUNDARY(Field3D, electrons["temperature"]);

  if (load_from_mesh){
    
  } else {
    phi = 3.0 * Te;
  }

  set(fields["phi"], phi);
  
}


void AdhocPotential::finally(const Options& state){
  AUTO_TRACE();

  const Options& allspecies = state["species"];

  
}

void AdhocPotential::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Tnorm = state["Tnorm"].as<BoutReal>();

  set_with_attrs(state["phi"], phi,
                 {{"time_dimension", "t"},
                  {"units", "V"},
                  {"conversion", Tnorm},
                  {"standard_name", "potential"},
                  {"long_name", "plasma potential"},
                  {"source", "adhoc_potential"}});
}


