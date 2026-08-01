#include "../include/sheath_boundary_fci.hxx"
#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/mesh.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/region.hxx>
#include <bout/sys/range.hxx>
#include <bout/utils.hxx>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using bout::globals::mesh;

namespace {


/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
///
BoutReal limitFree(BoutReal fm, BoutReal fc) {
  if (fm < fc) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }
  BoutReal fp = SQ(fc) / fm;
#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

} // namespace

BoutReal SheathBoundaryFci::ionSecondaryElectronEmissionGamma(BoutReal ion_energy) const {
  if (ion_ee_gamma_max < 0.0) {
    return 0.0;
  }
  if (ion_energy <= 0.0) {
    return 0.0;
  }

  return 0.5 * (1. + tanh((ion_energy - ion_ee_E_th) / (0.3 * ion_ee_E_th)))
         * ion_ee_gamma_max * std::pow(ion_energy / ion_ee_E_max, ion_ee_p)
         * std::exp(ion_ee_p * (1. - (ion_energy / ion_ee_E_max)));
}

SheathBoundaryFci::SheathBoundaryFci(std::string name, Options& alloptions, Solver*)
    : NamedComponent(name, {
                               readIfSet("species:{all_species}:charge"),
                               readIfSet("species:e:{e_whole_domain}"),
                               writeBoundary("species:e:{e_boundary}"),
                               readWrite("species:e:energy_source"),
                               writeBoundaryIfSet("species:e:{e_optional}"),
                               writeBoundaryReadInteriorIfSet("species:e:pressure"),
                               readIfSet("species:{ions}:adiabatic"),
                               readOnly("species:{ions}:AA"),
                               readWrite("species:{ions}:energy_source"),
                               writeBoundary("species:{ions}:{ion_boundary}"),
                               writeBoundaryReadInteriorIfSet("species:{ions}:pressure"),
                               writeBoundaryIfSet("species:{ions}:{ion_optional}"),
      }), yboundary(YBndryType::sheath, &alloptions, *mesh) {

  Options& options = alloptions[name];
  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];

  
  Ge = options["secondary_electron_coef"]
           .doc("Effective secondary electron emission coefficient")
           .withDefault(0.0);

  if ((Ge < 0.0) or (Ge > 1.0)) {
    throw BoutException("Secondary electron emission must be between 0 and 1 ({:e})", Ge);
  }

  ion_ee_gamma_max = options["ion_ee_gamma_max"]
                         .doc("Maximum ion induced secondary electron emission "
                              "coefficient. < 0 means off")
                         .withDefault(-1.0);
  ion_ee_E_th =
      options["ion_ee_E_th"]
          .doc("Energy threshold [eV] for ion induced secondary electron emission")
          .withDefault(50.0)
      / Tnorm;
  ion_ee_E_max =
      options["ion_ee_E_max"]
          .doc("Energy of maximum ion induced secondary electron emission [eV]")
          .withDefault(5e3)
      / Tnorm;

  ion_ee_p = options["ion_ee_p"]
                 .doc("Shape factor for ion induced secondary electron emission")
                 .withDefault(1.0);

  sin_alpha = options["sin_alpha"]
                  .doc("Sin of the angle between magnetic field line and wall surface. "
                       "Should be between 0 and 1")
                  .withDefault(Field3D(1.0));

  if (mesh->isFci()) {
    sin_alpha.applyBoundary();
    mesh->communicate(sin_alpha);
    sin_alpha.applyParallelBoundary("parallel_neumann_o1");
  }
  
  if ((min(sin_alpha) < 0.0) or (max(sin_alpha) > 1.0)) {
    throw BoutException("Range of sin_alpha must be between 0 and 1");
  }

  lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);

  diagnose = options["diagnose"].doc("Output diagnose variables?").withDefault<bool>(false);
  
  always_set_phi =
      options["always_set_phi"]
          .doc("Always set phi field? Default is to only modify if already set")
          .withDefault<bool>(false);

  // Read wall voltage, convert to normalised units
  wall_potential = options["wall_potential"]
                       .doc("Voltage of the wall [Volts]")
                       .withDefault(Field3D(0.0))
                   / Tnorm;
  // Convert to field aligned coordinates
  if (!mesh->isFci()) {
    wall_potential = toFieldAligned(wall_potential);
  } else {
    mesh->communicate(wall_potential);
    wall_potential.applyParallelBoundary("parallel_neumann_o1");
  }

  // Note: wall potential at the last cell before the boundary is used,
  // not the value at the boundary half-way between cells. This is due
  // to how twist-shift boundary conditions and non-aligned inputs are
  // treated; using the cell boundary gives incorrect results.

  floor_potential = options["floor_potential"]
                        .doc("Apply a floor to wall potential when calculating Ve?")
                        .withDefault<bool>(true);

  substitutePermissions("e_whole_domain", {"AA", "adiabatic"});
  substitutePermissions("e_boundary", {"density", "temperature"});
  substitutePermissions("e_optional", {"velocity", "momentum"});
  substitutePermissions("ion_boundary", {"density", "temperature"});
  substitutePermissions("ion_optional", {"velocity", "momentum"});
  setPermissions(always_set_phi ? writeBoundaryReadInteriorIfSet("fields:phi")
                                : writeBoundaryIfSet("fields:phi"));
}

void SheathBoundaryFci::transform_impl(GuardedOptions& state) {

  GuardedOptions allspecies = state["species"];
  GuardedOptions electrons = allspecies["e"];

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3D Ne = mesh->isFci()?
    floor(GET_NOBOUNDARY(Field3D, electrons["density"]).asField3DParallel(), 0.0)
    : toFieldAligned(floor(GET_NOBOUNDARY(Field3D, electrons["density"]).asField3DParallel(), 0.0));
  
  Field3D Te = mesh->isFci()?
    GET_NOBOUNDARY(Field3D, electrons["temperature"])
    :toFieldAligned(GET_NOBOUNDARY(Field3D, electrons["temperature"]));

  Field3D Pe;
  if (IS_SET_NOBOUNDARY(electrons["pressure"])) {
    Pe = mesh->isFci()?
      getNoBoundary<Field3D>(electrons["pressure"])
      :toFieldAligned(getNoBoundary<Field3D>(electrons["pressure"]));
  } else {
    Pe = Te * Ne.asField3DParallel();
  }
  
  // Ratio of specific heats
  const BoutReal electron_adiabatic =
      IS_SET(electrons["adiabatic"]) ? get<BoutReal>(electrons["adiabatic"]) : 5. / 3;

  // Mass, normalised to proton mass
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;

  // This is for applying boundary conditions
  Field3D Ve;
  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    Ve = mesh->isFci()?
      getNoBoundary<Field3D>(electrons["velocity"])
      : toFieldAligned(getNoBoundary<Field3D>(electrons["velocity"]));
  } else {
    Ve = zeroFrom(Ne);
  }

  Field3D NVe;
  if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
    NVe = mesh->isFci()?
      getNoBoundary<Field3D>(electrons["momentum"])
      : toFieldAligned(getNoBoundary<Field3D>(electrons["momentum"]));
  } else {
    NVe = zeroFrom(Ne);
  }
  

  // If Ve or NVe are not set, communicate in FCI. This is a work-around, as zeroFrom probably does not set the parallel slices (?)
  if (mesh->isFci()) {

    if (!IS_SET_NOBOUNDARY(electrons["velocity"])) {
      Ve.applyBoundary("neumann");
      mesh->communicate(Ve);
      Ve.applyParallelBoundary("parallel_neumann_o1");
    }

    if (!IS_SET_NOBOUNDARY(electrons["momentum"])) {
      NVe.applyBoundary("neumann");
      mesh->communicate(NVe);
      NVe.applyParallelBoundary("parallel_neumann_o1");
    }

    ASSERT2(Ne.hasParallelSlices());
    ASSERT2(Te.hasParallelSlices());
    ASSERT2(Pe.hasParallelSlices());
    ASSERT2(NVe.hasParallelSlices());
    ASSERT2(Ve.hasParallelSlices());


  }

  
  Coordinates* coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  // If phi is set, use free boundary condition
  // If phi not set, calculate assuming zero current
  Field3D phi;
  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    phi = phi.isFci()?
      getNoBoundary<Field3D>(state["fields"]["phi"])
      :toFieldAligned(getNoBoundary<Field3D>(state["fields"]["phi"]));
  } else {
    // Calculate potential phi assuming zero current
    // Note: This is equation (22) in Tskhakaya 2005, with I = 0

    struct IonBoundaryInfo {
      Field3D density;
      Field3D temperature;
      BoutReal mass;
      BoutReal charge;
      BoutReal adiabatic;
    };

    std::vector<IonBoundaryInfo> ion_species;

    // Need to sum  s_i Z_i C_i over all ion species
    //
    // To avoid looking up species for every grid point, this
    // loops over the boundaries once per species.
    Field3D ion_sum{zeroFrom(Ne)};

    if (mesh->isFci()) {
      mesh->communicate(ion_sum);
      ion_sum.applyParallelBoundary("parallel_neumann_o1");
    }
    
    phi = emptyFrom(Ne); // So phi is field aligned

    // Iterate through charged ion species
    for (auto& kv : allspecies.getChildren()) {
      GuardedOptions species = allspecies[kv.first];

      if ((kv.first == "e") or !IS_SET(species["charge"])
          or (get<BoutReal>(species["charge"]) == 0.0)) {
        continue; // Skip electrons and non-charged ions
      }

      const Field3D Ni = mesh->isFci()?
	floor(GET_NOBOUNDARY(Field3D, species["density"]).asField3DParallel(), 0.0)
	:toFieldAligned(floor(GET_NOBOUNDARY(Field3D, species["density"]), 0.0));
      const Field3D Ti = mesh->isFci()?
	GET_NOBOUNDARY(Field3D, species["temperature"])
	:toFieldAligned(GET_NOBOUNDARY(Field3D, species["temperature"]));
      const BoutReal Mi = GET_NOBOUNDARY(BoutReal, species["AA"]);
      const BoutReal Zi = GET_NOBOUNDARY(BoutReal, species["charge"]);

      const BoutReal adiabatic = IS_SET(species["adiabatic"])
                                     ? get<BoutReal>(species["adiabatic"])
                                     : 5. / 3; // Ratio of specific heats (ideal gas)

      if (ion_ee_gamma_max > 0.0) {
        ion_species.push_back({Ni, Ti, Mi, Zi, adiabatic});
      }

      

      ASSERT2(Ni.hasParallelSlices());
      ASSERT2(Ti.hasParallelSlices());
	

      yboundary.iter([&](auto& pnt) {	
	
	BoutReal s_i =
	  std::clamp(0.5 * (3. * pnt.current(Ni) / pnt.current(Ne) - pnt.prev(Ni) / pnt.prev(Ne)), 0.0, 1.0);

	if (!std::isfinite(s_i)) {
	  s_i = 1.0;
	}

	const BoutReal te = pnt.current(Te);
	const BoutReal ti = pnt.current(Ti);

	// Equation (9) in Tskhakaya 2005                                                                                                                                                                    

	BoutReal grad_ne = pnt.prev(Ne) - pnt.current(Ne);
	BoutReal grad_ni = pnt.prev(Ni) - pnt.current(Ni);
	  
	if (fabs(grad_ni) < 2e-3) {
	  grad_ni = grad_ne = 2e-3; // Remove kinetic correction term                                                                                                                                        
	}

	// Limit for e.g. Ni zero gradient                                                                                                                                                                   
	const BoutReal C_i_sq = std::clamp(
					   (adiabatic * ti + Zi * s_i * te * grad_ne / grad_ni) / Mi, 0., 100.);
	  
	// Note: Vzi = C_i * sin(α)                                                                                                                                                                          
	pnt.current(ion_sum) += s_i * Zi * pnt.prev(sin_alpha) * sqrt(C_i_sq);
	  
      }); 
    }
    
    
    phi.allocate();
    

    // ion_sum now contains  sum  s_i Z_i C_i over all ion species
    // at mesh->ystart and mesh->yend indices

    
    phi = 0.0;
    mesh->communicate(phi);
    ASSERT2(phi.hasParallelSlices());

    yboundary.iter([&](auto& pnt) {
      const auto& i = pnt.ind();

      if (pnt.current(Te) <= 0.0) {
	pnt.current(phi) = 0.0;
      } else {
	const BoutReal v_te = sqrt(pnt.current(Te) / (Me * TWOPI));
	const BoutReal prefactor = v_te * (1. - Ge);

	// Solve for sheath drop Δφ = φ - φ_wall, with optional ion-induced SEE:                                                                                                                             
	//   (1-Ge) v_te exp(-Δφ/Te) = ion_sum + Σ_i [γ_i(E_i(Δφ)) Γ_i/n_e]                                                                                                                                  
	BoutReal delta_phi = pnt.current(Te) * log(prefactor / pnt.current(ion_sum));

	if (ion_ee_gamma_max > 0.0) {
	  BoutReal delta_phi_guess = delta_phi;
	  constexpr int max_iter = 20;
	  constexpr BoutReal relax = 0.5;

	  for (int iter = 0; iter < max_iter; ++iter) {
	    BoutReal ion_emission_sum = 0.0;

	    for (const auto& ion : ion_species) {
	      const BoutReal te = pnt.current(Te);
	      const BoutReal ti = pnt.current(ion.temperature);
	      const BoutReal mi = ion.mass;
	      const BoutReal zi = ion.charge;
	      const BoutReal ad = ion.adiabatic;

	      BoutReal s_i = std::clamp(
					0.5 * (3. * pnt.current(ion.density) / pnt.current(Ne) - pnt.prev(ion.density) / pnt.prev(Ne)), 0.0,
					1.0);
		  
	      if (!std::isfinite(s_i)) {
		s_i = 1.0;
	      }

	      BoutReal grad_ne = pnt.prev(Ne) - pnt.current(Ne);
	      BoutReal grad_ni = pnt.prev(ion.density) - pnt.current(ion.density);

	      // Keep consistent with the ion_sum construction above
	      if (fabs(grad_ni) < 2e-3) {
		grad_ni = grad_ne = 2e-3;
	      }

	      const BoutReal C_i_sq = std::clamp(
						 (ad * ti + zi * s_i * te * grad_ne / grad_ni) / mi, 0., 100.);

	      const BoutReal ion_flux_over_ne =
		s_i * pnt.prev(sin_alpha) * sqrt(C_i_sq); // Γ_i / n_e

	      const BoutReal ion_energy =
		0.5 * (ti + mi * C_i_sq) + zi * delta_phi_guess;
	      ion_emission_sum +=
		ionSecondaryElectronEmissionGamma(ion_energy) * ion_flux_over_ne;
	    }
		
	    const BoutReal denom = ion_sum[i] + ion_emission_sum;
	    if (!(denom > 0.0) || !std::isfinite(denom)) {
	      break;
	    }

	    const BoutReal delta_phi_new = pnt.current(Te) * log(prefactor / denom);
	    if (!std::isfinite(delta_phi_new)) {
	      break;
	    }

	    const BoutReal err = fabs(delta_phi_new - delta_phi_guess);
	    const BoutReal scale = 1.0 + fabs(delta_phi_guess);
	    delta_phi_guess = (1. - relax) * delta_phi_guess + relax * delta_phi_new;

	    if (err < 1e-10 * scale) {
	      break;
	    }
	  }

	  delta_phi = delta_phi_guess;
	}

	pnt.current(phi) = delta_phi;
	    
      }

      const BoutReal phi_wall = pnt.current(wall_potential);
      pnt.current(phi) += phi_wall; // Add bias potential

      pnt.next(phi) = pnt.current(phi); // Constant into sheath
	  
    });
      
      
  }


  //////////////////////////////////////////////////////////////////
  // Electrons


  if (electrons.isSet("energy_source")) {
    electron_energy_source = mesh->isFci()?
      getNonFinal<Field3D>(electrons["energy_source"])
      :toFieldAligned(getNonFinal<Field3D>(electrons["energy_source"]));
  } else {
    electron_energy_source = zeroFrom(Ne);
  }


  yboundary.iter([&](auto& pnt) {
    const auto& i = pnt.ind();

    pnt.next(Ne) = limitFree(pnt.prev(Ne), pnt.current(Ne));
    pnt.next(Te) = limitFree(pnt.prev(Te), pnt.current(Te));
    //pnt.next(Pe) = limitFree(pnt.prev(Pe), pnt.current(Pe));
    pnt.next(Pe) = pnt.next(Ne) * pnt.next(Te);

    pnt.next(phi) = 2.0 * pnt.current(phi) - pnt.prev(phi);

    const BoutReal nesheath = 0.5 * (pnt.next(Ne) + pnt.current(Ne));
    const BoutReal tesheath = 0.5 * (pnt.next(Te) + pnt.current(Te)); // electron temperature
    const BoutReal phi_wall = pnt.current(wall_potential);

    const BoutReal phisheath =
      floor_potential ? floor(0.5 * (pnt.next(phi) + pnt.current(phi)),
			      phi_wall) // Electron saturation at phi = phi_wall                                                                                                                             
      : 0.5 * (pnt.current(phi) + pnt.next(phi));

    // Electron sheath heat transmission                                                                                                                                                                     
    const BoutReal gamma_e = floor(
				   (2 / (1. - Ge)) + ((phisheath - phi_wall) / floor(tesheath, 1e-5)), 0.0);
      
    // Electron velocity into sheath (< 0)                                                                                                                                                                   
    const BoutReal vesheath = (tesheath < 1e-10)
      ? 0.0
      : pnt.dir() * sqrt(tesheath / (TWOPI * Me)) * (1. - Ge)
      * exp(-(phisheath - phi_wall) / tesheath);

    pnt.next(Ve) = 2 * vesheath - pnt.prev(Ve);
    pnt.next(NVe) = 2. * Me * nesheath * vesheath - pnt.current(NVe);

    if (abs(pnt.offset()) == 1) { // Only subtract flux when the cell is actually in direct contact with the sheath
      
      BoutReal q = ((gamma_e - 1 - 1 / (electron_adiabatic - 1)) * tesheath
		    - 0.5 * Me * SQ(vesheath))
	* nesheath * vesheath;

      ASSERT2(q * pnt.dir() >= 0.0);
      
      // Multiply by cell area to get power                                                                                                                                                                    
      BoutReal flux;
      if (pnt.dir() > 0.0) {
	flux =  q * coord->cell_area_yhigh()[i];
      } else {
	flux =  q * coord->cell_area_ylow()[i];
      }
      
      
      // Divide by volume of cell to get energy loss rate (< 0)                                                                                                                                                
      const BoutReal power = flux / coord->cell_volume()[i];
      
#if CHECKLEVEL >= 1
      if (!std::isfinite(power)) {
	throw BoutException("Non-finite power at {} : Te {} Ne {} Ve {} phi {}, {}", i,
			    tesheath, nesheath, vesheath, pnt.current(phi), pnt.prev(phi));
      }
#endif
      
      electron_energy_source[i] -=  pnt.dir() * power;
    }
      
      
  });
  
  
  // Set electron density and temperature, now with boundary conditions
  // Note: Clear parallel slices because they do not contain correct boundary conditions
  setBoundary(electrons["density"], mesh->isFci()? Ne: fromFieldAligned(Ne));
  setBoundary(electrons["temperature"], mesh->isFci()? Te: fromFieldAligned(Te));
  setBoundary(electrons["pressure"], mesh->isFci()? Pe: fromFieldAligned(Pe));

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const

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
    Field3D Ni = mesh->isFci()?
      floor(getNoBoundary<Field3D>(species["density"]).asField3DParallel(), 0.0)
      :toFieldAligned(floor(getNoBoundary<Field3D>(species["density"]).asField3DParallel(), 0.0));
    Field3D Ti = mesh->isFci()?
      getNoBoundary<Field3D>(species["temperature"])
      :toFieldAligned(getNoBoundary<Field3D>(species["temperature"]));
    
    Field3D Pi = species.isSet("pressure")?
      (mesh->isFci()? getNoBoundary<Field3D>(species["pressure"]): toFieldAligned(getNoBoundary<Field3D>(species["pressure"])))
      : Ni * Ti.asField3DParallel();

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3D Vi = species.isSet("velocity")?
      (mesh->isFci()? getNoBoundary<Field3D>(species["velocity"]): toFieldAligned(getNoBoundary<Field3D>(species["velocity"])))
      : zeroFrom(Ni);
    Field3D NVi = species.isSet("momentum")?
      (mesh->isFci()? getNoBoundary<Field3D>(species["momentum"]): toFieldAligned(getNoBoundary<Field3D>(species["momentum"])))
      : Mi * Ni * Vi.asField3DParallel();

    if (mesh->isFci()) {
      
      if (!IS_SET_NOBOUNDARY(species["velocity"])) {
	Vi.applyBoundary("neumann");
	mesh->communicate(Vi);
	Vi.applyParallelBoundary("parallel_neumann_o1");
      }

      if (!IS_SET_NOBOUNDARY(species["momentum"])) {
	NVi.applyBoundary("neumann");
	mesh->communicate(NVi);
	NVi.applyParallelBoundary("parallel_neumann_o1");
      }

      ASSERT2(Ni.hasParallelSlices());
      ASSERT2(Ti.hasParallelSlices());
      ASSERT2(Pi.hasParallelSlices());
      ASSERT2(NVi.hasParallelSlices());
      ASSERT2(Vi.hasParallelSlices());
      
    }

    
    
    // Energy source will be modified in the domain
    Field3D energy_source =
        species.isSet("energy_source")?
      (mesh->isFci()? getNonFinal<Field3D>(species["energy_source"]): toFieldAligned(getNonFinal<Field3D>(species["energy_source"])))
            : zeroFrom(Ni);

  
    yboundary.iter([&](auto& pnt) {
      const auto& i = pnt.ind();
      
      pnt.next(Ni) = limitFree(pnt.prev(Ni), pnt.current(Ni));
      pnt.next(Ti) = limitFree(pnt.prev(Ti), pnt.current(Ti));
      //pnt.next(Pi) = limitFree(pnt.prev(Pi), pnt.current(Pi));	
      pnt.next(Pi) = pnt.next(Ni) * pnt.next(Ti);
      
      const BoutReal nesheath = floor(0.5 * (pnt.next(Ne) + pnt.current(Ne)), 1e-8);
      const BoutReal nisheath = floor(0.5 * (pnt.next(Ni) + pnt.current(Ni)), 1e-8);

      const BoutReal tesheath = floor(0.5 * (pnt.next(Te) + pnt.current(Te)), 1e-8); // electron temperature
      const BoutReal tisheath = floor(0.5 * (pnt.next(Ti) + pnt.current(Ti)), 1e-8); // electron temperature                                                                                                                 

      // Ion sheath heat transmission coefficient                                                                                                                                                            
      // Equation (22) in Tskhakaya 2005                                                                                                                                                                     
      // with                                                                                                                                                                                                
      //                                                                                                                                                                                                     
      // 1 / (1 + ∂_{ln n_e} ln s_i = s_i ∂_z n_e / ∂_z n_i                                                                                                                                                  
      // (from comparing C_i^2 in eq. 9 with eq. 20                                                                                                                                                          
      //                                                                                                                                                                                                     
      // Concentration                                                                                                                                                                                       
      const BoutReal s_i = std::clamp(nisheath / floor(nesheath, 1e-10), 0., 1.);
      BoutReal grad_ne = pnt.current(Ne) - nesheath;
      BoutReal grad_ni = pnt.current(Ni) - nisheath;
	
      if (fabs(grad_ni) < 1e-3) {
	grad_ni = grad_ne = 1e-3; // Remove kinetic correction term                                                                                                                                          
      }
	
      // Ion speed into sheath                                                                                                                                                                               
      // Equation (9) in Tskhakaya 2005                                                                                                                                                                      
      //                                                                                                                                                                                                     
      // Limit for e.g. Ni zero gradient                                                                                                                                                                     
      const BoutReal C_i_sq = std::clamp(
					 (adiabatic * tisheath + Zi * s_i * tesheath * grad_ne / grad_ni) / Mi, 0.,
					 100.);
	
      // Ion sheath heat transmission coefficient                                                                                                                                                            
      const BoutReal gamma_i = 2.5 + (0.5 * Mi * C_i_sq / tisheath);
	
      const BoutReal visheath = pnt.dir()>0?
	std::max(sqrt(C_i_sq), pnt.current(Vi))
	: std::min(-sqrt(C_i_sq), pnt.current(Vi));

      ASSERT2(visheath * pnt.dir() >= 0.0);
      
      // Set boundary conditions on flows                                                                                                                                                                    
      pnt.next(Vi) = 2. * visheath - pnt.current(Vi);
      pnt.next(NVi) = 2. * Mi * nisheath * visheath - pnt.current(NVi);
	
      if (abs(pnt.offset()) == 1) { // Only subtract flux when the cell is actually in direct contact with the sheath 
	
	BoutReal q =
	  floor(((gamma_i - 1 - 1 / (adiabatic - 1)) * tisheath - 0.5 * Mi * C_i_sq), 0.0)
	  * nisheath * visheath;

	ASSERT2(q * pnt.dir() >= 0.0);
	
	// Multiply by cell area to get power                                                                                                                                                                  
	BoutReal flux;
	if (pnt.dir() > 0) {
	  flux = q * coord->cell_area_yhigh()[i];
	} else {
	  flux = q * coord->cell_area_ylow()[i];
	}
	// Divide by volume of cell to get energy loss rate (< 0)                                                                                                                                              
	const BoutReal power = flux / coord->cell_volume()[i];
	ASSERT1(std::isfinite(power));
	ASSERT2(power * pnt.dir() >= 0.0);
	
	energy_source[i] -= pnt.dir() * power;
      }
	
      if (ion_ee_gamma_max > 0.0) {
	// Ion-induced secondary electron emission                                                                                                                                                           
	  
	const BoutReal phi_wall = pnt.current(wall_potential);
	const BoutReal phisheath = floor_potential
	  ? floor(0.5 * (pnt.next(phi) + pnt.current(phi)), phi_wall)
	  : 0.5 * (pnt.next(phi) + pnt.current(phi));
	  
	// Ion energy at the wall                                                                                                                                                                            
	const BoutReal Ei =
	  (0.5 * (tisheath + (Mi * SQ(visheath)))) + (Zi * (phisheath - phi_wall));
	  
	const BoutReal ion_ee_gamma = ionSecondaryElectronEmissionGamma(Ei);

	// Add flow of cold electrons into plasma                                                                                                                                                            
	BoutReal vesheath = 0.5 * (pnt.next(Ve) + pnt.current(Ve));
	const BoutReal nvesheath =
	  (nesheath * vesheath) - (ion_ee_gamma * nisheath * visheath);
	vesheath = nvesheath / nesheath;
	  
	pnt.next(Ve) = 2 * vesheath - pnt.current(Ve);
	pnt.next(NVe) = 2. * Me * nvesheath - pnt.current(NVe);
      }
	
	
	
    });
  
    setBoundary(species["density"], mesh->isFci()? Ni: fromFieldAligned(Ni));
    setBoundary(species["temperature"], mesh->isFci()? Ti: fromFieldAligned(Ti));
    setBoundary(species["pressure"], mesh->isFci()? Pi: fromFieldAligned(Pi));
    
    if (species.isSet("velocity")) {
      setBoundary(species["velocity"], mesh->isFci()? Vi: fromFieldAligned(Vi));
    }
    
    if (species.isSet("momentum")) {
      setBoundary(species["momentum"], mesh->isFci()? NVi: fromFieldAligned(NVi));
    }
  
    // Additional loss of energy through sheath
    // Note: Already includes previously set sources
    set(species["energy_source"], mesh->isFci()? energy_source: fromFieldAligned(energy_source));
    
  }
    
  // Add energy source (negative in cell next to sheath)
  // Note: already includes previously set sources
  set(electrons["energy_source"], mesh->isFci()? electron_energy_source: fromFieldAligned(electron_energy_source));
  
  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    setBoundary(electrons["velocity"], mesh->isFci()? Ve: fromFieldAligned(Ve));
  }
  if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
    setBoundary(electrons["momentum"], mesh->isFci()? NVe: fromFieldAligned(NVe));
  }
  if (always_set_phi or IS_SET_NOBOUNDARY(state["fields"]["phi"])) {    
    setBoundary(state["fields"]["phi"], mesh->isFci()? phi: fromFieldAligned(phi));
  }
}


void SheathBoundaryFci::outputVars(Options& state) {
  // Normalisations
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  
  
  if (diagnose) {
    set_with_attrs(
        state[std::string("electron_power")], electron_energy_source,
        {{"time_dimension", "t"},
         {"units", "Pa s^-1"},
         {"conversion", Tnorm * Nnorm * Omega_ci},
         {"source", "sheath_boundary_fci"}});
  }
  
}

