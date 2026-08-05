///////////////////////////////////////////////////////////////////////////
// LANGEVIN UPDATE - inside time loop
// Called each CFD timestep after turbulence->correct()
///////////////////////////////////////////////////////////////////////////

// Step 1: gradient of k - Y component (gravity direction)
gradK = fvc::grad(turbulence->k()());
dk_dy = gradK.component(vector::Y);

Info<< "dk_dy: min=" << min(dk_dy).value()
    << " max="        << max(dk_dy).value() << endl;

// Step 2: epsilon from k-epsilon solver
const volScalarField& epsilon_L =
    mesh.lookupObject<volScalarField>("epsilon");

Info<< "epsilon: min=" << min(epsilon_L).value()
    << " max="          << max(epsilon_L).value() << endl;

// Step 3: Lagrangian timescale T_L = 0.3*k/epsilon
T_L = dimensionedScalar("C_TL", dimless, C_TL)
    * turbulence->k()()
    / max
      (
          epsilon_L,
          dimensionedScalar("small",
          epsilon_L.dimensions(), small_k)
      );

Info<< "T_L: min=" << min(T_L).value()
    << " max="      << max(T_L).value() << endl;

// Step 4: drift correction - Y only
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

// Step 5: exponential decay
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

// Step 6: diffusion coefficient B
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

Info<< "B_coeff: min=" << min(B_coeff).value()
    << " max="          << max(B_coeff).value() << endl;

// Step 7: Langevin update loop
Random rndGen(label(runTime.timeIndex()));

forAll(mesh.cells(), cellI)
{
    scalar xi_x    = rndGen.GaussNormal<scalar>();
    scalar xi_y    = rndGen.GaussNormal<scalar>();
    scalar xi_z    = rndGen.GaussNormal<scalar>();
    scalar exp_i   = exp_dt[cellI];
    scalar B_i     = B_coeff[cellI];
    scalar drift_i = v_drift[cellI];

    // x: streamwise - no drift (periodic)
    up[cellI].x() = up[cellI].x() * exp_i
                  + B_i * xi_x;

    // y: wall-normal - drift applied (gravity direction)
    up[cellI].y() = up[cellI].y() * exp_i
                  + drift_i * (1.0 - exp_i)
                  + B_i * xi_y;

    // z: spanwise - no drift (periodic)
    up[cellI].z() = up[cellI].z() * exp_i
                  + B_i * xi_z;
}

// Step 8: BCs before u_seen
up.correctBoundaryConditions();

// Step 9: total velocity seen by particle
u_seen = U + up;

// Step 10: diagnostics
Info<< "Langevin diagnostics:" << nl
    << "  up_y   mean  = " << average(up.component(1)).value() << nl
    << "  up_y   std   = " << Foam::sqrt(average(
                                up.component(1)
                              * up.component(1))).value()       << nl
    << "  sqrt(2k/3)   = " << Foam::sqrt(average(
                                dimensionedScalar("c",
                                dimless, 2.0/3.0)
                              * turbulence->k()())).value()     << nl
    << "  u_seen max   = " << max(mag(u_seen)).value()         << nl
    << endl;

// Step 11: timestep check
scalar ratio_dt_TL =
    runTime.deltaT().value()
  / max(gMin(T_L), small_k);

if (ratio_dt_TL > 0.1)
{
    WarningInFunction
        << "dt/T_L = " << ratio_dt_TL
        << " exceeds 0.1. Reduce timestep for accurate Langevin."
        << endl;
}