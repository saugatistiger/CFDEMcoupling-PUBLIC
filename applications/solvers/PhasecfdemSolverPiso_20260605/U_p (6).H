///////////////////////////////////////////////////////////////////////////
// LANGEVIN UPDATE
// d50=0.01cm, dt=0.0005s, H_total=0.636cm (wall to free surface)
// Uses local k from simulation (physically correct)
// k clipped to k_inlet maximum (prevents corner spikes)
//
// References:
//   Thomson (1987) J. Fluid Mech. 180, 529-556  -- well-mixed condition
//   Pope (2000) Turbulent Flows, Cambridge       -- OU process, Appendix J
//   Van Rijn (1984) J. Hydraul. Eng. 110        -- reference height 2-5*D50
//   Rouse (1937)                                 -- concentration profile
///////////////////////////////////////////////////////////////////////////

const scalar dt_val   = runTime.deltaT().value();
const scalar d50_val  = D50.value();

///////////////////////////////////////////////////////////////////////////
// Smooth bed height — running average over 100 steps
// Prevents y_lang from oscillating each timestep which causes
// particles to switch abruptly between up=zero and full Langevin
// (this creates sharp discontinuities in the concentration profile)
///////////////////////////////////////////////////////////////////////////
static scalar y_bed_smooth = max(Bed_height_average, scalar(0.005));
y_bed_smooth = scalar(0.99) * y_bed_smooth
             + scalar(0.01) * max(Bed_height_average, scalar(0.005));
const scalar y_bed_ref = y_bed_smooth;

// Van Rijn (1984): reference height a = 2*D50 to 5*D50
// below this height: bed-load layer (damped Langevin)
// above this height: suspension zone (full Langevin)
// Use 2*D50 for low transport, 5*D50 for high transport
// Rouse number Z = ws/(kappa*u*): Z>1.2 -> 2D50, Z<0.8 -> 5D50
const scalar a_vanRijn_low  = scalar(0) * d50_val;  // low transport
const scalar a_vanRijn_high = scalar(5.0) * d50_val;  // high transport

// y_lang: Langevin activated above bed surface
const scalar y_lang   = y_bed_ref;
const scalar H_domain = max(scalar(63.6) * d50_val - y_bed_ref,
                            scalar(0.1));

///////////////////////////////////////////////////////////////////////////
// Physical scales from u_star
///////////////////////////////////////////////////////////////////////////
const scalar C_mu_val = 0.09;
const scalar sqrt_Cmu = Foam::sqrt(scalar(C_mu_val));
const scalar nu_val   = 0.01;  // kinematic viscosity [cm2/s]

scalar u_star_val = u_star.value();

Info<< "u_star (transportProperties) = " << u_star_val
    << "  Theta_x = " << Theta_average.x() << endl;

// k_inlet = u*^2 / sqrt(C_mu)
scalar k_inlet = u_star_val * u_star_val / sqrt_Cmu;

// Rouse number Z = ws/(kappa*u*) — determines reference height
scalar ws_val  = scalar(0.9);  // settling velocity [cm/s]
scalar Z_rouse = ws_val / (scalar(0.41) * max(u_star_val, scalar(1e-6)));

// Van Rijn (1984): blend reference height between 2*D50 and 5*D50
// using Rouse number: Z>1.2 -> 2D50 (near-bed), Z<0.8 -> 5D50 (full susp)
scalar blend    = max(scalar(0.0),
                  min(scalar(1.0),
                  (scalar(1.2) - Z_rouse) / scalar(0.4)));
scalar a_ref    = a_vanRijn_low
                + blend * (a_vanRijn_high - a_vanRijn_low);

Info<< "Rouse Z="     << Z_rouse
    << "  a_ref="     << a_ref
    << " cm (Van Rijn 1984: 2-5*D50)"
    << "  blend="     << blend  << endl;

///////////////////////////////////////////////////////////////////////////
// Step 1: Local k clipped to k_inlet
///////////////////////////////////////////////////////////////////////////
volScalarField k_local
(
    IOobject
    (
        "k_local",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    min
    (
        turbulence->k()(),
        dimensionedScalar("k_max",
        turbulence->k()().dimensions(), k_inlet)
    )
);

///////////////////////////////////////////////////////////////////////////
// Step 2: u_prime_field = sqrt(2*k_local/3)
///////////////////////////////////////////////////////////////////////////
volScalarField u_prime_field
(
    IOobject
    (
        "u_prime_field",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::sqrt
    (
        max
        (
            dimensionedScalar("c23", dimless, 2.0/3.0) * k_local,
            dimensionedScalar("small_k2", k_local.dimensions(), 1e-10)
        )
    )
);

scalar u_rms_inlet = Foam::sqrt(2.0*k_inlet/3.0);
scalar u_rms_avg   = average(u_prime_field).value();

Info<< "u_star_val  = " << u_star_val           << nl
    << "k_inlet     = " << k_inlet              << nl
    << "k_local max = " << max(k_local).value() << nl
    << "u_rms_avg   = " << u_rms_avg            << nl
    << "u_rms/ws    = " << u_rms_avg/ws_val     << endl;

///////////////////////////////////////////////////////////////////////////
// Step 3: epsilon clipped to physical range
// floor: small_eps with correct epsilon dimensions
// ceiling: log-law value at D50 from wall
///////////////////////////////////////////////////////////////////////////
const volScalarField& epsilon_raw =
    mesh.lookupObject<volScalarField>("epsilon");

scalar epsilon_max_val =
    Foam::pow(scalar(C_mu_val), scalar(0.75))
  * Foam::pow(k_inlet, scalar(1.5))
  / (scalar(0.41) * d50_val);

dimensionedScalar small_eps
(
    "small_eps",
    epsilon_raw.dimensions(),
    1e-10
);

volScalarField epsilon_L
(
    IOobject
    (
        "epsilon_L",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    min
    (
        max(epsilon_raw, small_eps),
        dimensionedScalar("eps_max",
        epsilon_raw.dimensions(), epsilon_max_val)
    )
);

Info<< "epsilon_L: min=" << min(epsilon_L).value()
    << " max="            << max(epsilon_L).value()
    << " eps_max="        << epsilon_max_val << endl;

///////////////////////////////////////////////////////////////////////////
// Step 4: T_L = C_TL * k_local / epsilon_L
// k_local used (not raw k)
// per-cell Kolmogorov floor: T_L >= 6*tau_eta = 6*sqrt(nu/epsilon)
// spatial smoothing of T_L reduces cell-boundary jumps in profile
///////////////////////////////////////////////////////////////////////////
T_L = dimensionedScalar("C_TL", dimless, C_TL)
    * k_local
    / max(epsilon_L, small_eps);

scalar T_L_max = H_domain / max(u_star_val, scalar(1e-6));
T_L_max        = max(T_L_max, scalar(10.0) * dt_val);

forAll(mesh.cells(), cellI)
{
    scalar tau_eta_i    = Foam::sqrt(
        nu_val / max(epsilon_L[cellI], scalar(1e-10))
    );
    scalar T_L_min_phys = scalar(6.0) * tau_eta_i;
    scalar T_L_min_num  = scalar(5.0) * dt_val;
    scalar T_L_min_i    = max(T_L_min_phys, T_L_min_num);
    T_L[cellI] = min(max(T_L[cellI], T_L_min_i), T_L_max);
}

// smooth T_L spatially — reduces cell-boundary jumps in concentration profile
T_L = fvc::average(fvc::interpolate(T_L));

// re-apply bounds after smoothing
forAll(mesh.cells(), cellI)
{
    scalar T_L_min_i = scalar(5.0) * dt_val;
    T_L[cellI] = min(max(T_L[cellI], T_L_min_i), T_L_max);
}
T_L.correctBoundaryConditions();

Info<< "T_L: min=" << min(T_L).value()
    << " max="      << max(T_L).value()
    << " T_L_max="  << T_L_max << endl;

///////////////////////////////////////////////////////////////////////////
// Step 5: B_coeff — horizontal diffusion coefficient
// exact OU variance: B = sqrt(C0*eps*T_L/2 * (1-exp(-2dt/T_L)))
///////////////////////////////////////////////////////////////////////////
dimensionedScalar small_TL
(
    "small_TL",
    T_L.dimensions(),
    1e-10
);

dimensionedScalar small_Bsq
(
    "small_Bsq",
    epsilon_L.dimensions() * T_L.dimensions(),
    1e-10
);

volScalarField B_coeff
(
    IOobject
    (
        "B_coeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::sqrt
    (
        max
        (
            dimensionedScalar("C0", dimless, C0)
            * epsilon_L * T_L
            / dimensionedScalar("two", dimless, 2.0)
            * (scalar(1.0) - Foam::exp
               (
                   -scalar(2.0) * runTime.deltaT()
                   / max(T_L, small_TL)
               )),
            small_Bsq
        )
    )
);

forAll(mesh.cells(), cellI)
{
    B_coeff[cellI] = min(B_coeff[cellI], u_prime_field[cellI]);
}

Info<< "B_coeff: min=" << min(B_coeff).value()
    << " max="          << max(B_coeff).value() << endl;

///////////////////////////////////////////////////////////////////////////
// Step 6: eps_t diagnostics at a_ref (Van Rijn reference height)
///////////////////////////////////////////////////////////////////////////
scalar eps_t_aref = scalar(0.41) * u_star_val * a_ref
                  * (scalar(1.0) - a_ref / H_domain);
scalar sigma_y_aref = Foam::sqrt(
    eps_t_aref / max(min(T_L).value(), scalar(1e-10))
);

Info<< "Van Rijn a_ref=" << a_ref
    << "  eps_t(a_ref)=" << eps_t_aref
    << "  sigma_y(a_ref)=" << sigma_y_aref
    << "  vs ws=" << ws_val << endl;

///////////////////////////////////////////////////////////////////////////
// Step 7: Sub-stepped Langevin update
// Zones (Van Rijn 1984):
//   below y_bed:         packed bed  — up = zero
//   y_bed to y_bed+a_ref: bed-load   — damped Langevin, no suspension drift
//   above y_bed+a_ref:   suspension  — full Langevin with Ito drift
//
// Smooth profile improvements:
//   - linear ramp of eps_t from 0 at bed to physical value at a_ref
//   - tanh soft clip instead of hard 3-sigma clip
//   - adaptive N_sub sub-steps for stability
///////////////////////////////////////////////////////////////////////////
scalar T_L_min_val = max(gMin(T_L), scalar(1e-10));
label N_sub = max(
    label(1),
    label(std::ceil(dt_val / (scalar(0.2) * T_L_min_val)))
);
N_sub = min(N_sub, label(20));
scalar dt_sub = dt_val / scalar(N_sub);

Info<< "Sub-stepping: N_sub=" << N_sub
    << "  dt_sub="   << dt_sub
    << "  dt/T_Lmin=" << dt_val/T_L_min_val
    << (dt_val/T_L_min_val < scalar(0.2) ? "  OK" : "  WARN>0.2")
    << endl;

label nClipped = 0;

forAll(mesh.cells(), cellI)
{
    scalar y_i       = mesh.C()[cellI].y();
    scalar u_prime_i = u_prime_field[cellI];
    scalar up_max_i  = scalar(3.0) * u_prime_i;

    if (y_i < y_lang)
    {
        // ── ZONE 1: packed bed — no motion ───────────────────────────
        up[cellI] = vector::zero;
    }
    else
    {
        scalar xi_i  = y_i - y_bed_ref;
        scalar T_L_i = T_L[cellI];
        scalar B_i   = B_coeff[cellI];

        if (xi_i <= a_ref)
        {
            // ── ZONE 2: bed-load layer (Van Rijn 1984: 0 to a_ref) ───
            // damped Langevin — small kicks, no suspension drift
            // linear ramp of B from 0 at bed to 0.3*B at a_ref
            // prevents abrupt transition at a_ref
            scalar frac   = xi_i / max(a_ref, scalar(1e-10));
            scalar B_bl   = frac * B_i * scalar(0.3);

            for (label isub = 0; isub < N_sub; isub++)
            {
                scalar xi_x   = rndGen.GaussNormal<scalar>();
                scalar xi_y   = rndGen.GaussNormal<scalar>();
                scalar xi_z   = rndGen.GaussNormal<scalar>();
                scalar exp_sub = Foam::exp(
                    -dt_sub / max(T_L_i, scalar(1e-10))
                );

                up[cellI].x() = up[cellI].x() * exp_sub + B_bl * xi_x;
                up[cellI].y() = up[cellI].y() * exp_sub + B_bl * xi_y;
                up[cellI].z() = up[cellI].z() * exp_sub + B_bl * xi_z;
            }

            // hard clip in bed-load — keep kicks small
            scalar up_bl_max = scalar(0.3) * u_prime_i;
            up[cellI].x() = min(max(up[cellI].x(),-up_bl_max), up_bl_max);
            up[cellI].y() = min(max(up[cellI].y(),-up_bl_max), up_bl_max);
            up[cellI].z() = min(max(up[cellI].z(),-up_bl_max), up_bl_max);
        }
        else
        {
            // ── ZONE 3: suspension zone (above a_ref) ────────────────
            // full Langevin with Ito drift (Thomson 1987)

            // parabolic eps_t with LINEAR RAMP near bed (0 to a_ref)
            // replaces abrupt floor at 1e-6 — gives smooth profile
            scalar eps_t_phys = scalar(0.41) * u_star_val * xi_i
                              * (scalar(1.0) - xi_i / H_domain);
            scalar ramp       = min(xi_i / max(a_ref, scalar(1e-10)),
                                    scalar(1.0));
            scalar eps_t_i    = max(eps_t_phys,
                                    ramp * eps_t_aref * scalar(0.1));

            // Ito drift: d(eps_t)/d(xi) — Thomson (1987) well-mixed
            // positive near bed (upward), zero at H/2, negative near surface
            scalar drift_y_i = scalar(0.41) * u_star_val
                             * (scalar(1.0) - scalar(2.0)*xi_i/H_domain);

            // 3-sigma wall-normal bound for tanh soft clip
            scalar up_y_sig_i = Foam::sqrt(
                max(eps_t_i / max(T_L_i, scalar(1e-10)), scalar(0.0))
            );
            scalar up_y_max_i = scalar(3.0) * up_y_sig_i;

            // sub-step loop
            for (label isub = 0; isub < N_sub; isub++)
            {
                scalar xi_x   = rndGen.GaussNormal<scalar>();
                scalar xi_y   = rndGen.GaussNormal<scalar>();
                scalar xi_z   = rndGen.GaussNormal<scalar>();

                scalar exp_sub = Foam::exp(
                    -dt_sub / max(T_L_i, scalar(1e-10))
                );

                // B for sub-step (horizontal)
                scalar B_sub = Foam::sqrt(
                    max(
                        scalar(C0) * epsilon_L[cellI] * T_L_i
                        / scalar(2.0)
                        * (scalar(1.0) - exp_sub * exp_sub),
                        scalar(0.0)
                    )
                );
                B_sub = min(B_sub, B_i);

                // B_y for sub-step (wall-normal, parabolic)
                scalar B_y_sub = Foam::sqrt(
                    max(
                        eps_t_i * (scalar(1.0) - exp_sub * exp_sub)
                        / max(T_L_i, scalar(1e-10)),
                        scalar(0.0)
                    )
                );
                B_y_sub = min(B_y_sub, B_sub);

                // x: streamwise
                up[cellI].x() = up[cellI].x() * exp_sub + B_sub * xi_x;
                up[cellI].x() = max(up[cellI].x(), -U[cellI].x());

                // y: wall-normal — Thomson (1987) Eq 30, no extra T_L
                up[cellI].y() = up[cellI].y() * exp_sub
                              + drift_y_i * (scalar(1.0) - exp_sub)
                              + B_y_sub * xi_y;

                // z: spanwise
                up[cellI].z() = up[cellI].z() * exp_sub + B_sub * xi_z;
            }
            // end sub-step loop

            // TANH SOFT CLIP — smoother than hard 3-sigma clip
            // saturates gradually instead of cutting off hard
            // produces smoother concentration profile near bed and surface
            if (up_y_max_i > scalar(1e-10))
            {
                up[cellI].x() = up_max_i
                    * Foam::tanh(up[cellI].x()
                    / max(up_max_i, scalar(1e-10)));
                up[cellI].y() = up_y_max_i
                    * Foam::tanh(up[cellI].y()
                    / max(up_y_max_i, scalar(1e-10)));
                up[cellI].z() = up_max_i
                    * Foam::tanh(up[cellI].z()
                    / max(up_max_i, scalar(1e-10)));
            }
        }
    }

    if (mag(up[cellI]) > up_max_i) { nClipped++; }
}

if (nClipped > 0)
{
    Info<< "up clipped: " << nClipped
        << "/" << mesh.nCells() << " cells"
        << " (" << 100.0*nClipped/mesh.nCells() << "%)" << endl;
}

///////////////////////////////////////////////////////////////////////////
// Step 8: boundary conditions
///////////////////////////////////////////////////////////////////////////
up.correctBoundaryConditions();

///////////////////////////////////////////////////////////////////////////
// Step 9: u_seen = U_mean + up
///////////////////////////////////////////////////////////////////////////
u_seen = U + up;

///////////////////////////////////////////////////////////////////////////
// NaN check
///////////////////////////////////////////////////////////////////////////
label nBad = 0;
forAll(mesh.cells(), cellI)
{
    if
    (
        !std::isfinite(u_seen[cellI].x()) ||
        !std::isfinite(u_seen[cellI].y()) ||
        !std::isfinite(u_seen[cellI].z())
    )
    {
        u_seen[cellI] = U[cellI];
        nBad++;
    }
}
if (nBad > 0)
{
    WarningInFunction
        << "u_seen: " << nBad
        << " non-finite values replaced with U_mean" << endl;
}

///////////////////////////////////////////////////////////////////////////
// Step 10: diagnostics
///////////////////////////////////////////////////////////////////////////
scalar up_y_mean = average(up.component(1)).value();
scalar up_y_std  = Foam::sqrt(
    average(up.component(1) * up.component(1)).value()
);

if (mag(up_y_mean) > scalar(0.1) * u_rms_avg)
{
    WarningInFunction
        << "up_y mean = " << up_y_mean
        << " not near zero — well-mixed condition violated." << nl
        << "Check parabolic eps_t profile and u_star value." << endl;
}

Info<< "Langevin diagnostics:" << nl
    << "  u_star_val  = " << u_star_val                              << nl
    << "  Z_rouse     = " << Z_rouse                                 << nl
    << "  a_ref       = " << a_ref << " cm (Van Rijn 2-5*D50)"      << nl
    << "  k_inlet     = " << k_inlet                                 << nl
    << "  u_rms_avg   = " << u_rms_avg                              << nl
    << "  N_sub       = " << N_sub                                   << nl
    << "  dt_sub      = " << dt_sub                                  << nl
    << "  T_L_min     = " << T_L_min_val                            << nl
    << "  up_y mean   = " << up_y_mean                              << nl
    << "  up_y std    = " << up_y_std                               << nl
    << "  up_y/u_rms  = " << up_y_std/max(u_rms_avg, scalar(1e-6)) << nl
    << "  up_y max    = " << max(up.component(1)).value()           << nl
    << "  up_y min    = " << min(up.component(1)).value()           << nl
    << "  u_seen max  = " << max(mag(u_seen)).value()               << nl
    << endl;

///////////////////////////////////////////////////////////////////////////
// Step 11: timestep stability check
///////////////////////////////////////////////////////////////////////////
scalar ratio = dt_sub / T_L_min_val;
if (ratio > scalar(0.2))
{
    WarningInFunction
        << "dt_sub/T_L = " << ratio
        << " exceeds 0.2 — increase N_sub or reduce dt" << nl
        << "dt_sub = "     << dt_sub
        << "  T_L_min = "  << T_L_min_val << endl;
}
