#include "../include/sheath_boundary_parallel.hxx"

#include <bout/output_bout_types.hxx>
#include <bout/yboundary_regions.hxx>

#include "bout/constants.hxx"
#include "bout/mesh.hxx"

using bout::globals::mesh;

namespace {
BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}

Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}
}


BoutReal smooth_step(BoutReal x, BoutReal f1, BoutReal f2) {

  if (x <= f1) {
    return 0.0;
  } else if (x >= f2) {
    return 1.0;
  } else {
    BoutReal t = (x-f1) / (f2 - f1);
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);    
  }
  
}

extern Options* tracking;
SheathBoundaryParallel::SheathBoundaryParallel(std::string name, Options &alloptions, Solver *)
  : yboundary(YBndryType::all, nullptr, *mesh){
  
  Options &options = alloptions[name];

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  
  Ge = options["secondary_electron_coef"]
           .doc("Effective secondary electron emission coefficient")
           .withDefault(0.0);

  if ((Ge < 0.0) or (Ge > 1.0)) {
    throw BoutException("Secondary electron emission must be between 0 and 1 ({:e})", Ge);
  }
  
  sin_alpha = options["sin_alpha"]
                  .doc("Sin of the angle between magnetic field line and wall surface. "
                       "Should be between 0 and 1")
                  .withDefault(1.0);


  bool dampen_low_density;
  BoutReal dampen_N_low,	dampen_N_high;

  dampen_low_density = options["dampen_low_density"]
          .doc("Dampen the sheath velocity to zero in low density regions?")
          .withDefault<bool>(false);

  dampen_N_low = options["dampen_N_low"].withDefault(3e17) / Nnorm;

  dampen_N_high = options["dampen_N_high"].withDefault(5e17) / Nnorm;


  

  

  if ((sin_alpha < 0.0) or (sin_alpha > 1.0)) {
    throw BoutException("Range of sin_alpha must be between 0 and 1");
  }

  always_set_phi =
      options["always_set_phi"]
          .doc("Always set phi field? Default is to only modify if already set")
          .withDefault<bool>(false);

  always_zero_current = options["always_zero_current"]
    .doc("Always set zero current?").withDefault<bool>(false);


  // Read wall voltage, convert to normalised units
  wall_potential = options["wall_potential"]
                       .doc("Voltage of the wall [Volts]")
                       .withDefault(Field3D(0.0))
                   / Tnorm;

  // Note: wall potential at the last cell before the boundary is used,
  // not the value at the boundary half-way between cells. This is due
  // to how twist-shift boundary conditions and non-aligned inputs are
  // treated; using the cell boundary gives incorrect results.

  sheath_extrapolate = options["sheath_extrapolate"]
                        .doc("Extrapolate values into the sheath? If not use neumann on all variables.")
                        .withDefault<bool>(true);
  
  floor_potential = options["floor_potential"]
                        .doc("Apply a floor to wall potential when calculating Ve?")
                        .withDefault<bool>(true);
}

void SheathBoundaryParallel::transform(Options &state) {

  Options& allspecies = state["species"];
  Options& electrons = allspecies["e"];

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3D Ne = toFieldAligned(floor(GET_NOBOUNDARY(Field3D, electrons["density"]), 0.0));
  Field3D Te = toFieldAligned(GET_NOBOUNDARY(Field3D, electrons["temperature"]));
  Field3D Pe = IS_SET_NOBOUNDARY(electrons["pressure"])
    ? toFieldAligned(getNoBoundary<Field3D>(electrons["pressure"]))
    : Te * Ne;

  // Ratio of specific heats
  const BoutReal electron_adiabatic =
      IS_SET(electrons["adiabatic"]) ? get<BoutReal>(electrons["adiabatic"]) : 5. / 3;

  // Mass, normalised to proton mass
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;

  // This is for applying boundary conditions
  Field3D Ve = IS_SET_NOBOUNDARY(electrons["velocity"])
    ? toFieldAligned(getNoBoundary<Field3D>(electrons["velocity"]))
    : zeroFrom(Ne);

  bool has_NVe = IS_SET_NOBOUNDARY(electrons["momentum"]);
  Field3D NVe;
  if (has_NVe) {
    NVe = toFieldAligned(getNoBoundary<Field3D>(electrons["momentum"]));
  }

  Coordinates *coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  // If phi is set, use free boundary condition
  // If phi not set, calculate assuming zero current
  Field3D phi;
  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    phi = toFieldAligned(getNoBoundary<Field3D>(state["fields"]["phi"]));
  } else {
    phi = zeroFrom(Ne); // So phi is field aligned
    always_zero_current = true; // Assume zero current to calculate phi on boundary
  }

  if (always_zero_current) {
    // Calculate potential phi assuming zero current
    // Note: This is equation (22) in Tskhakaya 2005, with I = 0

    // Need to sum  s_i Z_i C_i over all ion species
    //
    // To avoid looking up species for every grid point, this
    // loops over the boundaries once per species.
    Field3D ion_sum {zeroFrom(Ne)};

    // Iterate through charged ion species
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first];

      if ((kv.first == "e") or !IS_SET(species["charge"])
          or (get<BoutReal>(species["charge"]) == 0.0)) {
        continue; // Skip electrons and non-charged ions
      }

      const Field3D Ni = toFieldAligned(floor(GET_NOBOUNDARY(Field3D, species["density"]), 0.0));
      const Field3D Ti = toFieldAligned(GET_NOBOUNDARY(Field3D, species["temperature"]));
      const BoutReal Mi = GET_NOBOUNDARY(BoutReal, species["AA"]);
      const BoutReal Zi = GET_NOBOUNDARY(BoutReal, species["charge"]);

      const BoutReal adiabatic = IS_SET(species["adiabatic"])
                                     ? get<BoutReal>(species["adiabatic"])
                                     : 5. / 3; // Ratio of specific heats (ideal gas)

      yboundary.iter([&](auto& pnt) {
	
	const auto& i = pnt.ind();
	
	BoutReal s_i;
	if (sheath_extrapolate) {
	  BoutReal ratio = pnt.extrapolate_next_o2(Ni) / pnt.extrapolate_next_o2(Ni);
	  s_i = std::clamp(ratio, 1e-9, 1.0);
	} else {
	  s_i = pnt.current(Ni) / pnt.current(Ne);
	} 
	

	if (!std::isfinite(s_i)) {
	  s_i = 1.0;
	}
	BoutReal te = Te[i];
	BoutReal ti = Ti[i];
	
	// Equation (9) in Tskhakaya 2005
	BoutReal grad_ne;
	BoutReal grad_ni;
	
	if (sheath_extrapolate) {
	  grad_ne = pnt.extrapolate_grad_o2(Ne);
	  grad_ni = pnt.extrapolate_grad_o2(Ni);
	} else {
	  grad_ne = 1.0;
	  grad_ni = 1.0;
	}
	  
	// Note: Needed to get past initial conditions, perhaps
	// transients but this shouldn't happen in steady state
	if (fabs(grad_ni) < 1e-3) {
	  grad_ni = grad_ne = 1e-3; // Remove kinetic correction term
	}

	BoutReal C_i_sq =
	  clip((adiabatic * ti + Zi * s_i * te * grad_ne / grad_ni) / Mi, 0,
	       100); // Limit for e.g. Ni zero gradient
	
	// Note: Vzi = C_i * sin(α)
	BoutReal toadd = s_i * Zi * sin_alpha * sqrt(C_i_sq);
	if (legacy_match && pnt.dir() == 1) {
	  // sin_alpha missing
	  toadd = s_i * Zi * sqrt(C_i_sq);
	}
	ion_sum[i] += toadd;
      
      }); // end iter_regions
      
      
      phi = 0.0;
      
      // ion_sum now contains  sum  s_i Z_i C_i over all ion species
      // at mesh->ystart and mesh->yend indices
      yboundary.iter([&](auto& pnt) {      
	auto i = pnt.ind();     
	BoutReal thisphi;
	if (Te[i] <= 0.0) {
	  thisphi = 0.0;
	} else {
	  thisphi = Te[i] * log(sqrt(Te[i] / (Me * TWOPI)) * (1. - Ge) / ion_sum[i]);
	}
	
	thisphi += wall_potential[i];
	pnt.current(phi) = thisphi;
	pnt.next(phi) = thisphi;

      }); // end iter_regions
    }

  }
  //////////////////////////////////////////////////////////////////
  // Electrons

  Field3D electron_energy_source = electrons.isSet("energy_source")
    ? toFieldAligned(getNonFinal<Field3D>(electrons["energy_source"]))
    : zeroFrom(Ne);

  yboundary.iter([&](auto& pnt) {

    auto i = pnt.ind();

    // Free gradient of log electron density and temperature
    // Limited so that the values don't increase into the sheath
    // This ensures that the guard cell values remain positive
    // exp( 2*log(N[i]) - log(N[ip]) )
    if (sheath_extrapolate) {
      //pnt.limitFree(Ne);
      //pnt.limitFree(Te);
      //pnt.limitFree(Pe);
      pnt.extrapolate_next_o2(Ne);
      pnt.extrapolate_next_o2(Te);
      pnt.extrapolate_next_o2(Pe);
    } else {
      pnt.next(Ne) = pnt.current(Ne);
      pnt.next(Te) =	pnt.current(Te);
      pnt.next(Pe) =	pnt.current(Pe);
    }
      
      // Free boundary potential linearly extrapolated
    const BoutReal phiGradient = pnt.extrapolate_grad_o2(phi);
    pnt.neumann_o1(phi, phiGradient);
    
    const BoutReal nesheath = pnt.interpolate_boundary_o2(Ne);
    const BoutReal tesheath = pnt.interpolate_boundary_o2(Te);  // electron temperature
    const BoutReal phi_wall = pnt.current(wall_potential);
    
    const BoutReal phisheath = floor_potential ? floor(pnt.interpolate_boundary_o2(phi), phi_wall) // Electron saturation at phi = phi_wall
      : pnt.interpolate_boundary_o2(phi);

    // Electron sheath heat transmission
    const BoutReal gamma_e = floor(2 / (1. - Ge) + (phisheath - phi_wall) / floor(tesheath, 1e-5), 0.0);
    
    // Electron velocity into sheath (< 0)
    BoutReal vesheath = (tesheath < 1e-10) ?
      0.0 :
      pnt.dir() * sqrt(tesheath / (TWOPI * Me)) * (1. - Ge) * exp(-(phisheath - phi_wall) / tesheath);
    
    if (dampen_low_density) {
      vesheath = smooth_step(nesheath, dampen_N_low, dampen_N_high) * vesheath;
    }
    
      
      
    pnt.dirichlet_o2(Ve, vesheath);
    if (has_NVe) {
      pnt.dirichlet_o2(NVe, Me * nesheath * vesheath);
    }

    // Take into account the flow of energy due to fluid flow
    // This is additional energy flux through the sheath
    // Note: sign depends on sign of vesheath
    BoutReal q = ((gamma_e - 1 - 1 / (electron_adiabatic - 1)) * tesheath
		  - 0.5 * Me * SQ(vesheath))
      * nesheath * vesheath;
    
    // Multiply by cell area to get power
    BoutReal flux = 0.0;
    
    if (pnt.dir() < 0.0) {
      flux = q * coord->cell_area_ylow()[i];
    } else {
      flux = q * coord->cell_area_yhigh()[i];
    }

    // Divide by volume of cell to get energy loss rate (sign depending on vesheath)
    const BoutReal power = flux / coord->cell_volume()[i];
    
    electron_energy_source[i] -= pnt.dir() * power;
  
  }); // end iter_regions

  // Set electron density and temperature, now with boundary conditions
  setBoundary(electrons["density"], fromFieldAligned(Ne));
  setBoundary(electrons["temperature"], fromFieldAligned(Te));
  setBoundary(electrons["pressure"], fromFieldAligned(Pe));

  // Add energy source (negative in cell next to sheath)
  // Note: already includes previously set sources
  set(electrons["energy_source"], fromFieldAligned(electron_energy_source));

  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    setBoundary(electrons["velocity"], fromFieldAligned(Ve));
  }
  if (has_NVe) {
    setBoundary(electrons["momentum"], fromFieldAligned(NVe));
  }

  if (always_set_phi or IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    // Set the potential, including boundary conditions
    phi = fromFieldAligned(phi);
    //output.write("-> phi {}\n", phi(10, mesh->yend+1, 0));
    setBoundary(state["fields"]["phi"], phi);
  }

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Ion charge
    const BoutReal Zi =
        IS_SET(species["charge"]) ? get<BoutReal>(species["charge"]) : 0.0;

    if (Zi == 0.0) {
      continue; // Neutral -> skip
    }

    // Characteristics of this species
    const BoutReal Mi = get<BoutReal>(species["AA"]);

    const BoutReal adiabatic = IS_SET(species["adiabatic"])
                                   ? get<BoutReal>(species["adiabatic"])
                                   : 5. / 3; // Ratio of specific heats (ideal gas)

    // Density and temperature boundary conditions will be imposed (free)
    Field3D Ni = toFieldAligned(floor(getNoBoundary<Field3D>(species["density"]), 0.0));
    Field3D Ti = toFieldAligned(getNoBoundary<Field3D>(species["temperature"]));
    Field3D Pi = species.isSet("pressure")
      ? toFieldAligned(getNoBoundary<Field3D>(species["pressure"]))
      : Ni * Ti;

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3D Vi = species.isSet("velocity")
      ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
      : zeroFrom(Ni);
    Field3D NVi = species.isSet("momentum")
      ? toFieldAligned(getNoBoundary<Field3D>(species["momentum"]))
      : Mi * Ni * Vi;

    // Energy source will be modified in the domain
    Field3D energy_source = species.isSet("energy_source")
      ? toFieldAligned(getNonFinal<Field3D>(species["energy_source"]))
      : zeroFrom(Ni);

    yboundary.iter([&](auto& pnt) {


      auto i = pnt.ind();
      
      // Free gradient of log electron density and temperature
      // This ensures that the guard cell values remain positive
      // exp( 2*log(N[i]) - log(N[ip]) )
      if (sheath_extrapolate) {
	//pnt.limitFree(Ni);
	//pnt.limitFree(Ti);
	//pnt.limitFree(Pi);
	pnt.extrapolate_next_o2(Ni);
	pnt.extrapolate_next_o2(Ti);
	pnt.extrapolate_next_o2(Pi);
      } else {
	pnt.next(Ni) = pnt.current(Ni);
	pnt.next(Ti) = pnt.current(Ti);
	pnt.next(Pi) = pnt.current(Pi);
      }
      
      // Calculate sheath values at half-way points (cell edge)
      const BoutReal nesheath = pnt.interpolate_boundary_o2(Ne);
      const BoutReal nisheath = pnt.interpolate_boundary_o2(Ni);
      const BoutReal tesheath = floor(pnt.interpolate_boundary_o2(Te), 1e-5);  // electron temperature
      const BoutReal tisheath = floor(pnt.interpolate_boundary_o2(Ti), 1e-5);  // ion temperature
      
      // Ion sheath heat transmission coefficient
      // Equation (22) in Tskhakaya 2005
      // with 
      //
      // 1 / (1 + ∂_{ln n_e} ln s_i = s_i ∂_z n_e / ∂_z n_i
      // (from comparing C_i^2 in eq. 9 with eq. 20
      //
      //BoutReal s_i = (nesheath > 1e-5) ? nisheath / nesheath : 0.0; // Concentration ; upper_y
      BoutReal s_i = clip(nisheath / floor(nesheath, 1e-10), 0, 1); // Concentration ; lower_y
      if (legacy_match && pnt.dir() == -1){
	s_i = (nesheath > 1e-5) ? nisheath / nesheath : 0.0;
      }
      
      BoutReal grad_ne;
      BoutReal grad_ni;
      
      if (sheath_extrapolate) {
	grad_ne = pnt.extrapolate_grad_o2(Ne);
	grad_ni = pnt.extrapolate_grad_o2(Ni);
      } else {
	grad_ni = 1.0;
	grad_ne = 1.0;
      }
      
      if (fabs(grad_ni) < 1e-3) {
	grad_ni = grad_ne = 1e-3; // Remove kinetic correction term
      }
      
      // Ion speed into sheath
      // Equation (9) in Tskhakaya 2005
      //
      BoutReal C_i_sq =
	clip((adiabatic * tisheath + Zi * s_i * tesheath * grad_ne / grad_ni) / Mi,
	     0, 100); // Limit for e.g. Ni zero gradient
      
      if (dampen_low_density) {
	C_i_sq = smooth_step(nisheath, dampen_N_low, dampen_N_high) * C_i_sq;
      }
      
      const BoutReal visheath = pnt.dir() * sqrt(C_i_sq); // sign changes -> into sheath
      

	
      const BoutReal gamma_i = 2.5 + 0.5 * Mi * C_i_sq / tisheath; // + Δγ 


      // Set boundary conditions on flows
      pnt.dirichlet_o2(Vi, visheath);
      pnt.dirichlet_o2(NVi, Mi * nisheath * visheath);
      
      // Take into account the flow of energy due to fluid flow
      // This is additional energy flux through the sheath
      // Note: Sign depends on sign of visheath
      BoutReal q =
	((gamma_i - 1 - 1 / (adiabatic - 1)) * tisheath - 0.5 * C_i_sq * Mi)
	* nisheath * visheath;
      if (legacy_match and pnt.dir() == -1) {
	// Mi position switched with C_i_sq
	q =
	  ((gamma_i - 1 - 1 / (adiabatic - 1)) * tisheath - 0.5 * Mi * C_i_sq)
	  * nisheath * visheath;
      }
      
      if (q * pnt.dir() < 0.0) {
	q = 0.0;
      }
      
      // Multiply by cell area to get power
      BoutReal flux = 0.0;
      
      if (pnt.dir() < 0.0) {
	flux = q * coord->cell_area_ylow()[i];
      } else {
	flux = q * coord->cell_area_yhigh()[i];
	}
      
      
      // Divide by volume of cell to get energy loss rate (sign depending on vesheath)
      const BoutReal power = flux / coord->cell_volume()[i];
      
      ASSERT1(std::isfinite(power));
      ASSERT2(power * pnt.dir() >= 0.0);
      
      if (abs(pnt.offset()) == 1) {
	energy_source[pnt.ind()] -=
	  power * pnt.dir(); // Note: Sign negative because power * direction > 0
      }
    
    }); // end iter_regions

    // Finished boundary conditions for this species
    // Put the modified fields back into the state.
    setBoundary(species["density"], fromFieldAligned(Ni));
    setBoundary(species["temperature"], fromFieldAligned(Ti));
    setBoundary(species["pressure"], fromFieldAligned(Pi));

    if (species.isSet("velocity")) {
      setBoundary(species["velocity"], fromFieldAligned(Vi));
    }

    if (species.isSet("momentum")) {
      setBoundary(species["momentum"], fromFieldAligned(NVi));
    }
    if (tracking) {
      saveParallel(*tracking, fmt::format("NV{}_sheath", kv.first), NVi);
      saveParallel(*tracking, fmt::format("N{}_sheath",kv.first), Ni);
      saveParallel(*tracking, fmt::format("V{}_sheath", kv.first), Vi);
    }
    // Additional loss of energy through sheath
    // Note: Already includes previously set sources
    set(species["energy_source"], fromFieldAligned(energy_source));
  }

}
