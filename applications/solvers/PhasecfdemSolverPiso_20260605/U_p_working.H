///////////////////////////////////////////////////////////////////////////
// LANGEVIN UPDATE
// d50=0.01cm, dt=0.0001s, H=0.3cm
// k_inlet derived from prescribed u_star.value() (physically correct)
///////////////////////////////////////////////////////////////////////////

const scalar dt_val    = runTime.deltaT().value();
const scalar H_domain  = 0.003;    // H = 30*d50 in metres
                                    // 0.3cm = 0.003m

///////////////////////////////////////////////////////////////////////////
// k_inlet from prescribed u_star.value()
// k = u*^2 / sqrt(C_mu)
// this is the CORRECT turbulence level at inlet
// avoids using local k which can spike at corners
///////////////////////////////////////////////////////////////////////////
const scalar C_mu_val   = 0.09;
const scalar sqrt_Cmu   = Foam::sqrt(scalar(C_mu_val));   // = 0.3

scalar k_inlet          = u_star.value() * u_star.value() / sqrt_Cmu;
scalar u_prime_max      = (1.0/3.0) * Foam::sqrt(k_inlet);

Info<< "u_star.value()       = " << u_star.value()        << nl
    << "k_inlet      = " << k_inlet       << nl
    << "u_prime_max  = " << u_prime_max   << nl
    << "ws ~ 0.78 cm/s"                   << nl
    << "u'/ws        = " << u_prime_max/scalar(0.78) << endl;

///////////////////////////////////////////////////////////////////////////
// Step 1: gradient of k - Y component
// clip to physical range based on k_inlet
///////////////////////////////////////////////////////////////////////////
gradK = fvc::grad(turbulence->k()());
dk_dy = gradK.component(vector::Y);

scalar dk_dy_max = k_inlet / H_domain;   // physical scale k/H

forAll(mesh.cells(), cellI)
{
    dk_dy[cellI] = min(max(dk_dy[cellI], -dk_dy_max), dk_dy_max);
}
dk_dy.correctBoundaryConditions();

Info<< "dk_dy (clipped): min=" << min(dk_dy).value()
    << " max="                  << max(dk_dy).value()
    << " limit=+/-"             << dk_dy_max           << endl;

///////////////////////////////////////////////////////////////////////////
// Step 2: epsilon - clip extreme corner values
///////////////////////////////////////////////////////////////////////////
const volScalarField& epsilon_raw =
    mesh.lookupObject<volScalarField>("epsilon");

// epsilon_max based on k_inlet and minimum length scale
// epsilon_max ~ C_mu^0.75 * k_inlet^1.5 / (kappa * y_min)
// y_min ~ d50 = 0.01cm = 0.0001m
scalar y_min         = 0.0001;   // d50 in metres
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
    << " eps_max="        << epsilon_max_val          << endl;

///////////////////////////////////////////////////////////////////////////
// Step 3: T_L = 0.3*k/epsilon
// clip to [T_L_min, T_L_max]
// T_L_min = 5*dt  -> numerical stability
// T_L_max = H/u*  -> physical upper bound
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
scalar T_L_max = H_domain / max(u_star.value(), scalar(1e-6));

// safety: T_L_max >= 2*T_L_min
T_L_max = max(T_L_max, 2.0*T_L_min);

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
// clip v_drift to +/-sqrt(k_inlet) = u* scale
///////////////////////////////////////////////////////////////////////////
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
    (1.0/3.0) * T_L * dk_dy
);

scalar v_drift_max = Foam::sqrt(k_inlet);   // ~ u* scale

forAll(mesh.cells(), cellI)
{
    v_drift[cellI] = min(max(v_drift[cellI],
                             -v_drift_max), v_drift_max);
}

Info<< "v_drift (clipped): min=" << min(v_drift).value()
    << " max="                    << max(v_drift).value()
    << " limit=+/-"               << v_drift_max          << endl;

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
// capped at u_prime_max = (1/3)*sqrt(k_inlet)
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

// cap B at u_prime_max
forAll(mesh.cells(), cellI)
{
    B_coeff[cellI] = min(B_coeff[cellI], u_prime_max);
}

Info<< "B_coeff (clipped): min=" << min(B_coeff).value()
    << " max="                    << max(B_coeff).value()
    << " cap="                    << u_prime_max          << endl;

///////////////////////////////////////////////////////////////////////////
// Step 7: Langevin update
// clip each component to +/-u_prime_max
///////////////////////////////////////////////////////////////////////////
Random rndGen(label(runTime.timeIndex()));

label nClipped = 0;

forAll(mesh.cells(), cellI)
{
    scalar xi_x    = rndGen.GaussNormal<scalar>();
    scalar xi_y    = rndGen.GaussNormal<scalar>();
    scalar xi_z    = rndGen.GaussNormal<scalar>();
    scalar exp_i   = exp_dt[cellI];
    scalar B_i     = B_coeff[cellI];
    scalar drift_i = v_drift[cellI];

    // x: streamwise - memory + random
    up[cellI].x() = up[cellI].x() * exp_i
                  + B_i * xi_x;

    // y: wall-normal - memory + drift + random
    up[cellI].y() = up[cellI].y() * exp_i
                  + drift_i * (1.0 - exp_i)
                  + B_i * xi_y;

    // z: spanwise - memory + random
    up[cellI].z() = up[cellI].z() * exp_i
                  + B_i * xi_z;

    // clip each component to +/-u_prime_max
    up[cellI].x() = min(max(up[cellI].x(),
                            -u_prime_max), u_prime_max);
    up[cellI].y() = min(max(up[cellI].y(),
                            -u_prime_max), u_prime_max);
    up[cellI].z() = min(max(up[cellI].z(),
                            -u_prime_max), u_prime_max);

    if (mag(up[cellI]) > Foam::sqrt(scalar(3.0))*u_prime_max)
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
// FIX: check u_seen for NaN/Inf
// crash was in SS_turb_lift_mod receiving bad u_seen
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
        u_seen[cellI] = U[cellI];   // fallback to mean velocity
        nBad++;
    }
}

if (nBad > 0)
{
    WarningInFunction
        << "u_seen had " << nBad
        << " non-finite values -> replaced with U_mean"
        << endl;
}

///////////////////////////////////////////////////////////////////////////
// Step 10: diagnostics
///////////////////////////////////////////////////////////////////////////
Info<< "Langevin diagnostics:" << nl
    << "  u_star.value()       = " << u_star.value()                              << nl
    << "  k_inlet      = " << k_inlet                            << nl
    << "  u_prime_max  = " << u_prime_max                        << nl
    << "  up_y mean    = " << average(up.component(1)).value()   << nl
    << "  up_y std     = " << Foam::sqrt(average(
                                up.component(1)
                              * up.component(1))).value()         << nl
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
