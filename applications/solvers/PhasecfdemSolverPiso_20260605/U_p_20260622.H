///////////////////////////////////////////////////////////////////////////
// LANGEVIN UPDATE
// d50=0.01cm, dt=0.0001s, H=0.5cm
// Uses local k from simulation (physically correct)
// k clipped to k_inlet maximum (prevents corner spikes)
///////////////////////////////////////////////////////////////////////////

const scalar dt_val    = runTime.deltaT().value();
const scalar H_domain  = 0.5;    // flow depth: initial bed surface to free surface
const scalar y_bed_ref = 0.132;  // initial bed surface height from wall (top of lower mesh block)

///////////////////////////////////////////////////////////////////////////
// Physical scales from u_star
///////////////////////////////////////////////////////////////////////////
const scalar C_mu_val  = 0.09;
const scalar sqrt_Cmu  = Foam::sqrt(scalar(C_mu_val));

scalar u_star_val      = u_star.value();

// k_inlet = physical maximum k from inlet BC
// k = u*^2 / sqrt(C_mu)
scalar k_inlet         = u_star_val * u_star_val / sqrt_Cmu;

///////////////////////////////////////////////////////////////////////////
// Local k clipped to k_inlet
// uses local turbulence where physical
// prevents corner spikes
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
// u_prime_field = sqrt(2*k_local/3)
// PHYSICAL RMS turbulent fluctuation (not 1/3*sqrt(k))
// varies cell by cell using local k
// this is the correct turbulent intensity from theory:
//   u'_rms = sqrt(2k/3)
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
            dimensionedScalar("c23", dimless, 2.0/3.0)
          * k_local,
            dimensionedScalar("small",
            k_local.dimensions(), small_k)
        )
    )
);

// scalar reference values for diagnostics and limits
scalar u_rms_inlet     = Foam::sqrt(2.0*k_inlet/3.0);
scalar u_rms_avg       = average(u_prime_field).value();

Info<< "u_star_val    = " << u_star_val              << nl
    << "k_inlet       = " << k_inlet                 << nl
    << "k_local max   = " << max(k_local).value()    << nl
    << "k_local mean  = " << average(k_local).value()<< nl
    << "u_rms_inlet   = " << u_rms_inlet             << nl
    << "u_rms_avg     = " << u_rms_avg               << nl
    << "ws ~ 0.78 cm/s"                              << nl
    << "u_rms/ws      = " << u_rms_avg/0.78          << endl;

///////////////////////////////////////////////////////////////////////////
// Step 1: gradient of k - Y component
// clip using local k scale: dk/dy_max = k_local/H
///////////////////////////////////////////////////////////////////////////
gradK = fvc::grad(turbulence->k()());
dk_dy = gradK.component(vector::Y);

// physical gradient limit: k/H (integral scale)
forAll(mesh.cells(), cellI)
{
    scalar dk_dy_max_i = k_local[cellI] / H_domain;
    dk_dy[cellI] = min(max(dk_dy[cellI],
                           -dk_dy_max_i), dk_dy_max_i);
}
dk_dy.correctBoundaryConditions();

Info<< "dk_dy (clipped): min=" << min(dk_dy).value()
    << " max="                  << max(dk_dy).value() << endl;

///////////////////////////////////////////////////////////////////////////
// Step 2: epsilon - clip to physical range
// epsilon_max = C_mu^0.75 * k^1.5 / (kappa * y_min)
// y_min = d50 = 0.0001m
///////////////////////////////////////////////////////////////////////////
const volScalarField& epsilon_raw =
    mesh.lookupObject<volScalarField>("epsilon");

scalar y_min           = 0.01;
scalar epsilon_max_val =
    Foam::pow(scalar(C_mu_val), scalar(0.75))
  * Foam::pow(k_inlet, scalar(1.5))
  / (scalar(0.41) * y_min);

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
            dimensionedScalar("small",
            epsilon_raw.dimensions(), small_k)
        ),
        dimensionedScalar("eps_max",
        epsilon_raw.dimensions(), epsilon_max_val)
    )
);

Info<< "epsilon_L: min=" << min(epsilon_L).value()
    << " max="            << max(epsilon_L).value()
    << " eps_max="        << epsilon_max_val << endl;

///////////////////////////////////////////////////////////////////////////
// Step 3: T_L = 0.3*k/epsilon
// T_L_min = 5*dt  (numerical stability)
// T_L_max = H/u*  (physical upper bound)
///////////////////////////////////////////////////////////////////////////
T_L = dimensionedScalar("C_TL", dimless, C_TL)
    * turbulence->k()()
    / max
      (
          epsilon_L,
          dimensionedScalar("small",
          epsilon_L.dimensions(), small_k)
      );

scalar T_L_min = 5.0  * dt_val;
scalar T_L_max = H_domain / max(u_star_val, scalar(1e-6));
T_L_max        = max(T_L_max, 2.0*T_L_min);

forAll(mesh.cells(), cellI)
{
    T_L[cellI] = min(max(T_L[cellI], T_L_min), T_L_max);
}
T_L.correctBoundaryConditions();

Info<< "T_L: min="    << min(T_L).value()
    << " max="         << max(T_L).value()
    << " T_L_min="     << T_L_min
    << " T_L_max="     << T_L_max
    << " exp range=("  << Foam::exp(-dt_val/T_L_max)
    << ","             << Foam::exp(-dt_val/T_L_min) << ")" << endl;

///////////////////////////////////////////////////////////////////////////
// Step 4: drift correction - Y only
// v_drift = drift_coeff * T_L * dk/dy
// drift_coeff = 0.05 (reduced from 1/3 to fix upward bias)
// clip to local u_prime (not global scalar)
///////////////////////////////////////////////////////////////////////////
scalar drift_coeff = 1.0/3.0;

volScalarField v_drift
(
    IOobject
    (
        "v_drift",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    drift_coeff * T_L * dk_dy
);

// clip using local u_prime
forAll(mesh.cells(), cellI)
{
    scalar v_drift_max_i = u_prime_field[cellI];
    v_drift[cellI] = min(max(v_drift[cellI],
                             -v_drift_max_i), v_drift_max_i);
}

// remove global mean to prevent net upward bias when dk/dy > 0 throughout domain
scalar v_drift_mean_val = average(v_drift).value();
forAll(mesh.cells(), cellI)
{
    v_drift[cellI] -= v_drift_mean_val;
}

Info<< "v_drift mean removed: " << v_drift_mean_val              << nl
    << "v_drift (zero-mean): min=" << min(v_drift).value()
    << " max="                     << max(v_drift).value() << endl;

///////////////////////////////////////////////////////////////////////////
// Step 5: exponential decay
///////////////////////////////////////////////////////////////////////////
volScalarField exp_dt
(
    IOobject
    (
        "exp_dt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::exp
    (
        -runTime.deltaT()
        / max
          (
              T_L,
              dimensionedScalar("small",
              T_L.dimensions(), small_k)
          )
    )
);

Info<< "exp_dt: min=" << min(exp_dt).value()
    << " max="         << max(exp_dt).value() << endl;

///////////////////////////////////////////////////////////////////////////
// Step 6: diffusion coefficient B
// B naturally gives std(up) = sqrt(2k/3) = u_rms
// cap at local u_prime (one RMS per step maximum)
///////////////////////////////////////////////////////////////////////////
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
                   / max
                     (
                         T_L,
                         dimensionedScalar("small",
                         T_L.dimensions(), small_k)
                     )
               )),
            dimensionedScalar("small",
            epsilon_L.dimensions()*T_L.dimensions(), small_k)
        )
    )
);

// cap B at local u_prime (physically motivated)
forAll(mesh.cells(), cellI)
{
    B_coeff[cellI] = min(B_coeff[cellI],
                         u_prime_field[cellI]);
}

Info<< "B_coeff: min=" << min(B_coeff).value()
    << " max="          << max(B_coeff).value() << endl;

///////////////////////////////////////////////////////////////////////////
// Step 7: Langevin update
// clip using LOCAL u_prime (not global scalar)
// 3-sigma physical limit per cell
///////////////////////////////////////////////////////////////////////////
Random rndGen(label(runTime.timeIndex()));

// near-bed Lagrangian diffusivity floor
// Eulerian eps_t -> 0 at wall, but Lagrangian particle diffusivity does not:
// eddies from y+ > 30 sweep down and entrain particles at the bed.
// Without a floor, B_y -> 0 near bed traps particles and creates a mound.
// Floor = mixing-length value evaluated at y+ = 30 (log-layer start).
scalar nu_val      = turbulence->nu()().internalField()[0];
scalar y_inner     = 30.0 * nu_val / max(u_star_val, scalar(1e-6));
y_inner            = min(y_inner, H_domain * scalar(0.5));
scalar eps_t_inner = scalar(0.41) * u_star_val * y_inner
                   * (scalar(1.0) - y_inner / H_domain);

Info<< "eps_t_inner (y+=30 floor): " << eps_t_inner
    << " at y_inner=" << y_inner << endl;

label nClipped = 0;

forAll(mesh.cells(), cellI)
{
    scalar xi_x      = rndGen.GaussNormal<scalar>();
    scalar xi_y      = rndGen.GaussNormal<scalar>();
    scalar xi_z      = rndGen.GaussNormal<scalar>();
    scalar exp_i     = exp_dt[cellI];
    scalar B_i       = B_coeff[cellI];
    scalar drift_i   = v_drift[cellI];

    // local u_prime for this cell
    scalar u_prime_i = u_prime_field[cellI];

    // 3-sigma limit (99.7% of Gaussian distribution)
    scalar up_max_i  = 3.0 * u_prime_i;

    // wall-normal B for y-component: match parabolic mixing-length diffusivity
    // eps_t(xi) = kappa * u_star * xi * (1 - xi/H)  [open channel]
    // xi = height above initial bed surface (not above wall)
    // floor at y+=30 prevents B_y -> 0 trapping particles at the bed
    scalar y_i         = mesh.C()[cellI].y();
    scalar xi_i        = max(y_i - y_bed_ref, scalar(0.0));
    scalar eps_t_i     = max(
        scalar(0.41) * u_star_val * xi_i
        * (scalar(1.0) - xi_i / H_domain),
        eps_t_inner
    );
    scalar B_y_i   = Foam::sqrt(
        max(
            eps_t_i * (scalar(1.0) - exp_i * exp_i)
            / max(T_L[cellI], scalar(1e-10)),
            scalar(0.0)
        )
    );
    B_y_i = min(B_y_i, B_i);  // cannot exceed isotropic B

    // x: streamwise - memory + random (isotropic)
    up[cellI].x() = up[cellI].x() * exp_i
                  + B_i * xi_x;

    // y: wall-normal - memory + drift + mixing-length-damped random
    up[cellI].y() = up[cellI].y() * exp_i
                  + drift_i * (1.0 - exp_i)
                  + B_y_i * xi_y;

    // z: spanwise - memory + random (isotropic)
    up[cellI].z() = up[cellI].z() * exp_i
                  + B_i * xi_z;

    // clip to 3*u_prime_local (physically motivated)
    scalar up_y_max_i = 3.0 * Foam::sqrt(max(eps_t_i
                      / max(T_L[cellI], scalar(1e-10)), scalar(0.0)));
    up[cellI].x() = min(max(up[cellI].x(),
                            -up_max_i),   up_max_i);
    up[cellI].y() = min(max(up[cellI].y(),
                            -up_y_max_i), up_y_max_i);
    up[cellI].z() = min(max(up[cellI].z(),
                            -up_max_i),   up_max_i);

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

// Step 8: BCs
up.correctBoundaryConditions();

// Step 9: u_seen = U_mean + up
u_seen = U + up;

///////////////////////////////////////////////////////////////////////////
// NaN check - prevents crash in SS_turb_lift_mod
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
        << " non-finite -> replaced with U_mean"
        << endl;
}

///////////////////////////////////////////////////////////////////////////
// Step 10: diagnostics
// well-mixed: up_y mean ~ 0
// correct intensity: up_y std ~ u_rms_avg
///////////////////////////////////////////////////////////////////////////
scalar up_y_mean = average(up.component(1)).value();
scalar up_y_std  = Foam::sqrt
                   (
                       average
                       (
                           up.component(1)
                         * up.component(1)
                       ).value()
                   );

// well-mixed check
if (mag(up_y_mean) > 0.1*u_rms_avg)
{
    WarningInFunction
        << "up_y mean = " << up_y_mean
        << " not near zero! Well-mixed violated." << nl
        << "Consider reducing drift_coeff from "
        << drift_coeff << endl;
}

Info<< "Langevin diagnostics:" << nl
    << "  u_star_val   = " << u_star_val                          << nl
    << "  k_inlet      = " << k_inlet                            << nl
    << "  u_rms_inlet  = " << u_rms_inlet                        << nl
    << "  u_rms_avg    = " << u_rms_avg                          << nl
    << "  drift_coeff  = " << drift_coeff                        << nl
    << "  up_y mean    = " << up_y_mean                          << nl
    << "  up_y std     = " << up_y_std                           << nl
    << "  up_y/u_rms   = " << up_y_std/max(u_rms_avg,scalar(1e-6))<< nl
    << "  up_y max     = " << max(up.component(1)).value()        << nl
    << "  up_y min     = " << min(up.component(1)).value()        << nl
    << "  up     max   = " << max(mag(up)).value()                << nl
    << "  u_seen max   = " << max(mag(u_seen)).value()            << nl
    << endl;

///////////////////////////////////////////////////////////////////////////
// Step 11: timestep check
///////////////////////////////////////////////////////////////////////////
scalar T_L_min_val = max(gMin(T_L), small_k);
scalar ratio       = dt_val / T_L_min_val;

if (ratio > 0.2)
{
    WarningInFunction
        << "dt/T_L = "   << ratio
        << " exceeds 0.2" << nl
        << "dt = "        << dt_val
        << " T_L_min = "  << T_L_min_val
        << endl;
}
