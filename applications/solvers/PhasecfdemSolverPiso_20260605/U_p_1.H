///////////////////////////////////////////////////////////////////////////
// LANGEVIN UPDATE
// d50=0.01cm, dt=0.0005s, H_total=0.636cm (wall to free surface)
// Uses local k from simulation (physically correct)
// k clipped to k_inlet maximum (prevents corner spikes)
//
// References:
//   Thomson (1987) J. Fluid Mech. 180, 529-556  -- well-mixed condition
//   Pope (2000) Turbulent Flows, Cambridge       -- OU process, Appendix J
//   Rouse (1937)                                 -- concentration profile
///////////////////////////////////////////////////////////////////////////

const scalar dt_val    = runTime.deltaT().value();

// y_bed_ref: current bed surface from deltaY (no artificial shift)
const scalar y_bed_ref = max(Bed_height_average, scalar(0.005));

// y_lang: Langevin applied to all cells above bed surface
const scalar y_lang    = y_bed_ref;
const scalar H_domain  = max(scalar(63.6*D50.value()) - y_bed_ref,
                             scalar(0.1));

///////////////////////////////////////////////////////////////////////////
// Physical scales from u_star
///////////////////////////////////////////////////////////////////////////
const scalar C_mu_val  = 0.09;
const scalar sqrt_Cmu  = Foam::sqrt(scalar(C_mu_val));
const scalar nu_val    = 0.01;  // kinematic viscosity [cm2/s]

scalar u_star_val = u_star.value();

// Dynamic u_star override DISABLED:
// tau_average is computed from velocity gradients inside the dense
// particle bed where particle drag suppresses fluid velocity,
// giving u_star ~5x too low.
// Use u_star from transportProperties (correctly set per-case).
Info<< "u_star (transportProperties) = " << u_star_val
    << "  Theta_x = " << Theta_average.x() << endl;

// k_inlet = physical maximum k from inlet BC: k = u*^2 / sqrt(C_mu)
scalar k_inlet = u_star_val * u_star_val / sqrt_Cmu;

///////////////////////////////////////////////////////////////////////////
// Step 1: Local k clipped to k_inlet
// Prevents corner spikes from contaminating Langevin fluctuations
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
// Physical RMS turbulent fluctuation — varies cell by cell
// From isotropic turbulence: u'_rms = sqrt(2k/3)
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
            dimensionedScalar("small_k2",
            k_local.dimensions(), 1e-10)
        )
    )
);

scalar u_rms_inlet = Foam::sqrt(2.0*k_inlet/3.0);
scalar u_rms_avg   = average(u_prime_field).value();

Info<< "u_star_val    = " << u_star_val              << nl
    << "k_inlet       = " << k_inlet                 << nl
    << "k_local max   = " << max(k_local).value()    << nl
    << "k_local mean  = " << average(k_local).value()<< nl
    << "u_rms_inlet   = " << u_rms_inlet             << nl
    << "u_rms_avg     = " << u_rms_avg               << nl
    << "ws ~ 0.9 cm/s"                               << nl
    << "u_rms/ws      = " << u_rms_avg/0.9           << endl;

///////////////////////////////////////////////////////////////////////////
// Step 3: epsilon - clip to physical range
// epsilon_max = C_mu^0.75 * k^1.5 / (kappa * y_min)  [log-law at D50]
// FIX: use small_eps with correct epsilon dimensions (not small_k)
///////////////////////////////////////////////////////////////////////////
const volScalarField& epsilon_raw =
    mesh.lookupObject<volScalarField>("epsilon");

scalar y_min           = D50.value();  // reference height = D50
scalar epsilon_max_val =
    Foam::pow(scalar(C_mu_val), scalar(0.75))
  * Foam::pow(k_inlet, scalar(1.5))
  / (scalar(0.41) * y_min);

// FIX: correct epsilon dimensions for floor value
dimensionedScalar small_eps
(
    "small_eps",
    epsilon_raw.dimensions(),  // cm2/s3 — correct for epsilon
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
        max
        (
            epsilon_raw,
            small_eps              // FIX: correct dimensions and value
        ),
        dimensionedScalar("eps_max",
        epsilon_raw.dimensions(), epsilon_max_val)
    )
);

Info<< "epsilon_L: min=" << min(epsilon_L).value()
    << " max="            << max(epsilon_L).value()
    << " eps_max="        << epsilon_max_val << endl;

///////////////////////////////////////////////////////////////////////////
// Step 4: T_L = C_TL * k_local / epsilon_L
// FIX 1: use k_local (clipped) not raw turbulence->k()()
// FIX 2: use small_eps (correct dimensions) not small_k
// FIX 3: per-cell Kolmogorov floor (physical lower bound)
// T_L_max = H/u* (integral scale upper bound)
///////////////////////////////////////////////////////////////////////////

// FIX: use k_local and small_eps
T_L = dimensionedScalar("C_TL", dimless, C_TL)
    * k_local
    / max(epsilon_L, small_eps);

scalar T_L_max = H_domain / max(u_star_val, scalar(1e-6));
T_L_max        = max(T_L_max, scalar(10.0) * dt_val);

// FIX: per-cell Kolmogorov floor — physical lower bound
// Inertial subrange requires T_L >> tau_eta = sqrt(nu/epsilon)
forAll(mesh.cells(), cellI)
{
    // Kolmogorov timescale at this cell
    scalar tau_eta_i    = Foam::sqrt(
        nu_val / max(epsilon_L[cellI], scalar(1e-10))
    );

    // physical minimum: 6*tau_eta
    scalar T_L_min_phys = 6.0 * tau_eta_i;

    // numerical minimum: 5*dt for OU stability (dt/T_L < 0.2)
    scalar T_L_min_num  = 5.0 * dt_val;

    // take larger of physical and numerical floors
    scalar T_L_min_i    = max(T_L_min_phys, T_L_min_num);

    T_L[cellI] = min(max(T_L[cellI], T_L_min_i), T_L_max);
}
T_L.correctBoundaryConditions();

Info<< "T_L: min="    << min(T_L).value()
    << " max="         << max(T_L).value()
    << " T_L_max="     << T_L_max
    << " exp range=("  << Foam::exp(-dt_val/T_L_max)
    << ","             << Foam::exp(-dt_val/min(T_L).value())
    << ")" << endl;

///////////////////////////////////////////////////////////////////////////
// Step 5: B_coeff — diffusion coefficient (horizontal directions)
// From exact OU variance (Pope 2000, Appendix J, Eq J.41):
//   B = sqrt( C0*epsilon*T_L/2 * (1 - exp(-2*dt/T_L)) )
// FIX: correct dimension objects throughout
///////////////////////////////////////////////////////////////////////////

// FIX: correct small dimension for T_L (time dimensions)
dimensionedScalar small_TL
(
    "small_TL",
    T_L.dimensions(),
    1e-10
);

// FIX: correct small dimension for B^2 = epsilon*T_L
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
                   / max(T_L, small_TL)  // FIX: correct dimensions
               )),
            small_Bsq                    // FIX: correct dimensions
        )
    )
);

// cap B at local u_prime (physically motivated)
forAll(mesh.cells(), cellI)
{
    B_coeff[cellI] = min(B_coeff[cellI], u_prime_field[cellI]);
}

Info<< "B_coeff: min=" << min(B_coeff).value()
    << " max="          << max(B_coeff).value() << endl;

///////////////////////////////////////////////////////////////////////////
// Step 6: parabolic eps_t diagnostics at D50 and 2*D50
///////////////////////////////////////////////////////////////////////////
scalar xi_D50     = scalar(D50.value());
scalar eps_t_D50  = scalar(0.41) * u_star_val * xi_D50
                  * (scalar(1.0) - xi_D50 / H_domain);
scalar xi_2D50    = scalar(2.0 * D50.value());
scalar eps_t_2D50 = scalar(0.41) * u_star_val * xi_2D50
                  * (scalar(1.0) - xi_2D50 / H_domain);

Info<< "y_bed_ref = " << y_bed_ref
    << "  H_domain = " << H_domain                          << nl
    << "eps_t(1*D50)=" << eps_t_D50
    << "  eps_t(2*D50)=" << eps_t_2D50                     << nl
    << "sigma_y(2*D50)=" << Foam::sqrt(eps_t_2D50
       / max(min(T_L).value(), scalar(1e-10)))
    << "  vs ws=0.9"                                        << endl;

///////////////////////////////////////////////////////////////////////////
// Step 7: Sub-stepped Langevin update
// FIX: N_sub sub-steps per CFD timestep
//      ensures dt_sub/T_L <= 0.2 near bed
// FIX: adaptive N_sub based on minimum T_L
// FIX: RNG — use rndGen declared in solver (not re-seeded here)
// FIX: drift term WITHOUT extra T_L factor (Thomson 1987, Eq 30)
// FIX: clipping AFTER all sub-steps complete
///////////////////////////////////////////////////////////////////////////

// adaptive N_sub: ceil(dt / (0.2 * T_L_min))
scalar T_L_min_val = max(gMin(T_L), scalar(1e-10));
label N_sub = max(
    label(1),
    label(std::ceil(dt_val / (scalar(0.2) * T_L_min_val)))
);
N_sub = min(N_sub, label(20));  // cap at 20 sub-steps
scalar dt_sub = dt_val / scalar(N_sub);

Info<< "Langevin sub-stepping:"
    << "  N_sub="      << N_sub
    << "  dt_sub="     << dt_sub
    << "  dt/T_L_min=" << dt_val/T_L_min_val
    << (dt_val/T_L_min_val < 0.2 ? "  OK" : "  WARNING>0.2")
    << endl;

label nClipped = 0;

forAll(mesh.cells(), cellI)
{
    scalar y_i       = mesh.C()[cellI].y();
    scalar u_prime_i = u_prime_field[cellI];
    scalar up_max_i  = 3.0 * u_prime_i;

    if (y_i < y_lang)
    {
        // ── BELOW BED SURFACE: packed bed, no Langevin ───────────────
        up[cellI] = vector::zero;
    }
    else
    {
        // ── SUSPENSION ZONE: full Langevin above bed ─────────────────
        scalar xi_i  = y_i - y_bed_ref;
        scalar T_L_i = T_L[cellI];
        scalar B_i   = B_coeff[cellI];

        // parabolic eps_t at this height — constant across sub-steps
        // eps_t -> 0 at bed and free surface, peaks at H/2
        scalar eps_t_i = max(
            scalar(0.41) * u_star_val * xi_i
            * (scalar(1.0) - xi_i / H_domain),
            scalar(1e-6)
        );

        // Ito drift: d(eps_t)/d(xi) — Thomson (1987) well-mixed condition
        // for parabolic eps_t = kappa*u*xi(1-xi/H):
        //   d(eps_t)/d(xi) = kappa*u*(1 - 2*xi/H)
        // FIX: correct form WITHOUT extra T_L factor
        // positive near bed (upward), zero at H/2, negative near surface
        scalar drift_y_i = scalar(0.41) * u_star_val
                         * (scalar(1.0) - scalar(2.0) * xi_i / H_domain);

        // 3-sigma wall-normal bound (computed once, applied after sub-loop)
        scalar up_y_max_i = 3.0 * Foam::sqrt(
            max(eps_t_i / max(T_L_i, scalar(1e-10)),
                scalar(0.0))
        );

        // ── sub-step loop ─────────────────────────────────────────────
        for (label isub = 0; isub < N_sub; isub++)
        {
            // fresh independent Gaussian random numbers each sub-step
            scalar xi_x = rndGen.GaussNormal<scalar>();
            scalar xi_y = rndGen.GaussNormal<scalar>();
            scalar xi_z = rndGen.GaussNormal<scalar>();

            // exponential decay for sub-step dt_sub (exact OU)
            scalar exp_sub = Foam::exp(
                -dt_sub / max(T_L_i, scalar(1e-10))
            );

            // B for sub-step: exact OU variance (Pope 2000, App J)
            scalar B_sub = Foam::sqrt(
                max(
                    scalar(C0) * epsilon_L[cellI] * T_L_i
                    / scalar(2.0)
                    * (scalar(1.0) - exp_sub * exp_sub),
                    scalar(0.0)
                )
            );
            B_sub = min(B_sub, B_i);

            // B_y for sub-step: parabolic eps_t variance
            scalar B_y_sub = Foam::sqrt(
                max(
                    eps_t_i * (scalar(1.0) - exp_sub * exp_sub)
                    / max(T_L_i, scalar(1e-10)),
                    scalar(0.0)
                )
            );
            B_y_sub = min(B_y_sub, B_sub);

            // x: streamwise Langevin update
            up[cellI].x() = up[cellI].x() * exp_sub
                          + B_sub * xi_x;

            // prevent back-flow exceeding mean velocity
            up[cellI].x() = max(up[cellI].x(), -U[cellI].x());

            // y: wall-normal — Thomson (1987) Eq 30
            // drift_y_i*(1-exp_sub): correct form, no extra T_L
            up[cellI].y() = up[cellI].y() * exp_sub
                          + drift_y_i * (scalar(1.0) - exp_sub)
                          + B_y_sub * xi_y;

            // z: spanwise Langevin update
            up[cellI].z() = up[cellI].z() * exp_sub
                          + B_sub * xi_z;
        }
        // ── end sub-step loop ─────────────────────────────────────────

        // FIX: clip to 3-sigma AFTER all sub-steps complete (not inside loop)
        up[cellI].x() = min(max(up[cellI].x(), -up_max_i),   up_max_i);
        up[cellI].y() = min(max(up[cellI].y(), -up_y_max_i), up_y_max_i);
        up[cellI].z() = min(max(up[cellI].z(), -up_max_i),   up_max_i);
    }

    if (mag(up[cellI]) > up_max_i)
    {
        nClipped++;
    }
}

if (nClipped > 0)
{
    Info<< "up clipped: " << nClipped
        << "/" << mesh.nCells() << " cells"
        << " (" << 100.0*nClipped/mesh.nCells() << "%)"
        << endl;
}

///////////////////////////////////////////////////////////////////////////
// Step 8: boundary conditions
///////////////////////////////////////////////////////////////////////////
up.correctBoundaryConditions();

///////////////////////////////////////////////////////////////////////////
// Step 9: u_seen = U_mean + up
// This is what all force models use (DiFeliceDrag, gradPForce,
// viscForce, SS_turb_lift) — must be set to u_seen in couplingProperties
///////////////////////////////////////////////////////////////////////////
u_seen = U + up;

///////////////////////////////////////////////////////////////////////////
// NaN check — prevents crash in SS_turb_lift_mod
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
        << " non-finite values replaced with U_mean"
        << endl;
}

///////////////////////////////////////////////////////////////////////////
// Step 10: diagnostics
// well-mixed check: <up_y> ~ 0
// intensity check:  std(up_y) ~ u_rms_avg
///////////////////////////////////////////////////////////////////////////
scalar up_y_mean = average(up.component(1)).value();
scalar up_y_std  = Foam::sqrt
                   (
                       average
                       (
                           up.component(1) * up.component(1)
                       ).value()
                   );

// well-mixed condition check (Thomson 1987)
if (mag(up_y_mean) > 0.1 * u_rms_avg)
{
    WarningInFunction
        << "up_y mean = " << up_y_mean
        << " not near zero — well-mixed condition violated." << nl
        << "Check parabolic eps_t profile and u_star value." << endl;
}

Info<< "Langevin diagnostics:" << nl
    << "  u_star_val   = " << u_star_val                              << nl
    << "  k_inlet      = " << k_inlet                                 << nl
    << "  u_rms_inlet  = " << u_rms_inlet                            << nl
    << "  u_rms_avg    = " << u_rms_avg                              << nl
    << "  N_sub        = " << N_sub                                   << nl
    << "  dt_sub       = " << dt_sub                                  << nl
    << "  T_L_min      = " << T_L_min_val                            << nl
    << "  up_y mean    = " << up_y_mean                              << nl
    << "  up_y std     = " << up_y_std                               << nl
    << "  up_y/u_rms   = " << up_y_std/max(u_rms_avg, scalar(1e-6)) << nl
    << "  up_y max     = " << max(up.component(1)).value()           << nl
    << "  up_y min     = " << min(up.component(1)).value()           << nl
    << "  up     max   = " << max(mag(up)).value()                   << nl
    << "  u_seen max   = " << max(mag(u_seen)).value()               << nl
    << endl;

///////////////////////////////////////////////////////////////////////////
// Step 11: timestep stability check
///////////////////////////////////////////////////////////////////////////
scalar ratio = dt_sub / T_L_min_val;

if (ratio > 0.2)
{
    WarningInFunction
        << "dt_sub/T_L = " << ratio
        << " exceeds 0.2 — increase N_sub or reduce dt" << nl
        << "dt_sub = "     << dt_sub
        << "  T_L_min = "  << T_L_min_val
        << endl;
}
