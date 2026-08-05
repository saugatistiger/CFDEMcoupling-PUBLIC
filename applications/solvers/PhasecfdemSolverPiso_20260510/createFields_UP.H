///////////////////////////////////////////////////////////////////////////
// Langevin fields - INITIALIZATION ONLY
// Include this in createFields.H
///////////////////////////////////////////////////////////////////////////

Info<< "Initializing Langevin turbulent dispersion fields" << endl;

// Turbulent fluctuation velocity
// READ_IF_PRESENT: preserves memory term on restart
volVectorField up
(
    IOobject
    (
        "up",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
);

// Wall-normal gradient of k (Y = gravity direction)
volScalarField dk_dy
(
    IOobject
    (
        "dk_dy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 1, -2, 0, 0), 0.0)
);

// gradient of k field
volVectorField gradK
(
    IOobject
    (
        "gradK",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector(0,0,0))
);

// Lagrangian timescale
volScalarField T_L
(
    IOobject
    (
        "T_L",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0)
);

// Total velocity seen by particle = U_mean + up
volVectorField u_seen
(
    IOobject
    (
        "u_seen",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
);

// Constants
const scalar C0      = 2.1;
const scalar C_TL    = 0.3;
const scalar small_k = 1e-10;