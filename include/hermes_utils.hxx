#pragma once
#ifndef HERMES_UTILS_H
#define HERMES_UTILS_H

#include "bout/bout_enum_class.hxx"
#include <algorithm>

#include "bout/traits.hxx"

inline BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}


inline BoutReal softFloor(BoutReal value, BoutReal min) {
  BoutReal limvalue = std::max(value, 0.0);
  return limvalue + min * exp(-limvalue / min);
}



inline Field3D softFloor(const Field3D& var, BoutReal f, const std::string& rgn = "RGN_NOY") {
  checkData(var);
  if (var.hasParallelSlices()) {
    Field3DParallel result;
    result.allocate();

    BOUT_FOR(i, var.getRegion(rgn)) {
      const auto iyp = i.yp();
      const auto iym = i.ym();
      result[i] = softFloor(var[i], f);
      result.yup()[iyp] = softFloor(var.yup()[iyp], f);
      result.yup()[iym] = softFloor(var.ydown()[iym], f);      
    }    
    return result.asField3D();
    
  } else {
    Field3D result;
    result.allocate();
    BOUT_FOR(d, var.getRegion(rgn)) { result[d] = softFloor(var[d], f); }
    return result;
  }
}




template<typename T, typename = bout::utils::EnableIfField<T>>
inline T clamp(const T& var, BoutReal lo, BoutReal hi, const std::string& rgn = "RGN_ALL") {
  checkData(var);
  T result = copy(var);

  BOUT_FOR(d, var.getRegion(rgn)) {
    if (result[d] < lo) {
      result[d] = lo;
    } else if (result[d] > hi) {
      result[d] = hi;
    }
  }

  return result;
}

/// Enum that identifies the type of a species: electron, ion, neutral
BOUT_ENUM_CLASS(SpeciesType, electron, ion, neutral);

/// Identify species name string as electron, ion or neutral
inline SpeciesType identifySpeciesType(const std::string& species) {
  if (species == "e") {
    return SpeciesType::electron;
  } else if ((species == "i") or
             species.find(std::string("+")) != std::string::npos) {
    return SpeciesType::ion;
  }
  // Not electron or ion -> neutral
  return SpeciesType::neutral;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
Ind3D indexAt(const T& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

#endif // HERMES_UTILS_H
