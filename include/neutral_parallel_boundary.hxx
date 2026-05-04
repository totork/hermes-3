#pragma once
#ifndef NEUTRAL_PARALLEL_BOUNDARY_H
#define NEUTRAL_PARALLEL_BOUNDARY_H

#include "component.hxx"

#include <bout/yboundary_regions.hxx>

struct NeutralParallelBoundary : public Component {


  NeutralParallelBoundary(std::string name, Options& alloptions, Solver*);

  void transform(Options &state) override;

private:
  BoutReal boundary_value;
  YBoundary yboundary;

  Field3D fromFieldAligned(const Field3D& f) {
    if (f.isFci()) {
      return f;
    }
    return ::fromFieldAligned(f);
  }
  Field3D toFieldAligned(const Field3D& f) {
    if (f.isFci()) {
      return f;
    }
    return ::toFieldAligned(f);
  }
  template <class F>
  void iter_regions(const F& f) {
    yboundary.iter_regions(f);
  }
  
};

namespace {
RegisterComponent<NeutralParallelBoundary> registercomponentneutralparallelboundary("neutral_parallel_boundary");
}

#endif
