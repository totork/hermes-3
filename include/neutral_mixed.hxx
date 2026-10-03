
#pragma once
#ifndef NEUTRAL_MIXED_H
#define NEUTRAL_MIXED_H

#include <memory>
#include <string>

#include <bout/invert_laplace.hxx>

#include "component.hxx"
#include <bout/yboundary_regions.hxx>
#include <bout/derivs.hxx>



/// Evolve density, parallel momentum and pressure
/// for a neutral gas species with cross-field diffusion
struct NeutralMixed : public Component {
  ///
  /// @param name     The name of the species e.g. "h"
  /// @param options  Top-level options. Settings will be taken from options[name]
  /// @param solver   Time-integration solver to be used
  NeutralMixed(const std::string& name, Options& options, Solver *solver);
  
  /// Modify the given simulation state
  void transform(Options &state) override;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

  /// Preconditioner
  void precon(const Options &state, BoutReal gamma) override;
private:
  std::string name;  ///< Species name
  YBoundary yboundary;
  std::shared_ptr<FCI::dagp_fv> dagp;
  Field3D Nn, Pn, NVn; // Density, pressure and parallel momentum
  Field3D Vn; ///< Neutral parallel velocity
  Field3D Tn; ///< Neutral temperature
  Field3D Nnlim, Pnlim, logPnlim, Vnlim, Tnlim; // Limited in regions of low density
  bool isMMS;
  Field3D Pn_solver;
  bool use_eos;
  Field3D ones;
  bool viscous_heating;
  BoutReal AA; ///< Atomic mass (proton = 1)
  BoutReal n_lowsource, T_lowsource, lowsource_scale;
  bool lowsource_balance;
  std::string balance_species;
  BoutReal low_N_lim;
  Field3D Dnn; ///< Diffusion coefficient
  Field3D DnnNn, DnnPn, DnnNVn;
  bool disable_Dnn;
  BoutReal temperature_floor;
  bool sheath_ydown, sheath_yup;
  bool dissipative;
  bool upwinding;
  BoutReal density_floor; ///< Minimum Nn used when dividing NVn by Nn to get Vn.
  BoutReal pressure_floor; ///< Minimum Pn used when dividing Pn by Nn to get Tn.
  bool exponential_source;
  BoutReal flux_limit; ///< Diffusive flux limit
  bool LP_limit;
  BoutReal LP_speed;
  Field3D lambdaLP;
  BoutReal limit_length;
  BoutReal diffusion_limit;    ///< Maximum diffusion coefficient
  bool inherited_T;
  BoutReal include_cond;
  Field3D anomalous_conduction;
  bool disable_ddt;
  bool parallel_dirichlet;
  bool neutral_viscosity; ///< include viscosity?
  bool neutral_conduction; ///< Include heat conduction?
  bool evolve_momentum; ///< Evolve parallel momentum?
  bool evolve_pressure;
  bool freeze_low_density;
  BoutReal freeze_density_value;
  bool use_finite_difference;
  Field3D initial_Tn;
  Field3D kappa_n, eta_n; ///< Neutral conduction and viscosity
  BoutReal neutral_lmax;
  int diffusion_mode;
  bool precondition {true}; ///< Enable preconditioner?
  bool lax_flux; ///< Use Lax flux for advection terms
  std::unique_ptr<Laplacian> inv; ///< Laplacian inversion used for preconditioning
  bool simplified_diffusion;
  BoutReal Dnn_last_update, Dnn_update_every;
  Field3D Dnn_cached, kappa_n_cached, eta_n_cached;
  BoutReal cond_factor;
  Field3D density_source, pressure_source; ///< External input source
  Field3D Sn, Sp, Snv; ///< Particle, pressure and momentum source
  Field3D sound_speed; ///< Sound speed for use with Lax flux
  BoutReal sound_speed_Tfloor;
  bool output_ddt; ///< Save time derivatives?
  bool diagnose; ///< Save additional diagnostics?
  bool output_transport;
  Field3D Nh_up, Nh_down;
  int precon_mode;
  // Flow diagnostics
  Field3D pf_adv_perp_xlow, pf_adv_perp_ylow, pf_adv_par_ylow;
  Field3D mf_adv_perp_xlow, mf_adv_perp_ylow, mf_adv_par_ylow;
  Field3D mf_visc_perp_xlow, mf_visc_perp_ylow, mf_visc_par_ylow;
  Field3D ef_adv_perp_xlow, ef_adv_perp_ylow, ef_adv_par_ylow;
  Field3D ef_cond_perp_xlow, ef_cond_perp_ylow, ef_cond_par_ylow;

  std::string equi_species;
  bool core_equilibriate;
  BoutReal tau_eq, n_thresh, delta_n, t_thresh, delta_t;
  BoutReal lowsource_width;

  const Field3D smooth_low_sourceterm(const Field3D& f, const BoutReal lowvalue,
                              const BoutReal scalefactor,
                              const BoutReal width) {
    Field3D result = 0.0;

    ASSERT1(width > 0.0);
    
    BOUT_FOR(i, f.getRegion("RGN_NOY")) {
      const BoutReal diff = f[i] - lowvalue;
      const BoutReal s = -diff / width; // = (lowvalue - f) / width
      
      if (s <= -1.0) {
	// Comfortably above threshold: no source
	result[i] = 0.0;
      } else if (s >= 1.0) {
	// Comfortably below threshold: original linear ramp
	result[i] = -diff / scalefactor;
      } else {
	// Smooth C2 transition region: quintic Hermite blend
	const BoutReal x = 0.5 * (s + 1.0);
	const BoutReal x3 = x * x * x;
	const BoutReal x4 = x3 * x;
	const BoutReal p = 2.0 * x3 - x4; // 0 at x=0, 1 at x=1, matching derivatives
	
	result[i] = (width / scalefactor) * p;
      }
    }
    
    return result;
  }
  
  
  const Field3D Grad_x(Field3D& a) {
    Mesh* mesh = a.getMesh();
    Coordinates* coord = mesh->getCoordinates();
    return DDX(a) / sqrt(coord->g_11);
  }


  const Field3D Grad_z(Field3D& a) {
    Mesh* mesh = a.getMesh();
    Coordinates* coord = mesh->getCoordinates();
    return DDZ(a) / sqrt(coord->g_33);
  }
  
  Field3D Div_a_Grad_perp_neutrals(Field3D a, Field3D b, std::shared_ptr<FCI::dagp_fv> dagp_op, bool upwind, int mode) {

    if (mode == 1) {
      return (*dagp_op)(a,b,upwind);
    } else if (mode ==2 ) {
      return Div_a_Grad_perp_curv(a, b);
    }
    
  }
};

namespace {
RegisterComponent<NeutralMixed> registersolverneutralmixed("neutral_mixed");
}

#endif // NEUTRAL_MIXED_H
