#pragma once
#ifndef ADHOC_POTENTIAL_H
#define ADHOC_POTENTIAL_H

#include "component.hxx"

struct AdhocPotential : public Component {
  AdhocPotential(std::string name, Options &options, Solver*UNUSED(solver));

  void transform(Options& state) override;

  void finally(Options& state) override;

  void outputVars(Options& state) override;

private:
  
  Field3D phi;
  Field3D Te;
  
  bool load_from_mesh;
  

};

namespace {
  RegisterComponent<AdhocPotential> registercomponentrelaxpotential("adhoc_potential");
}

#endif 
