    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//===============================
// particle interaction modelling
//===============================

    Info<< "\nReading momentum exchange field Ksl\n" << endl;
    volScalarField Ksl
    (
        IOobject
        (
            "Ksl",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
    );

    Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "\nCreating density field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
        //dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );

    Info<< "Reading particle velocity field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//===============================
 Info<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
 volScalarField gradAlpha
 (
     IOobject
     (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     mag(fvc::grad(voidfraction))
 );

//#   include "createPhi.H"
#ifndef createPhi_H
#define createPhi_H
Info<< "Equating Alpha = voidfraction" << endl;
volScalarField alpha
 (
     IOobject
     (
        "alpha",
        runTime.timeName(),
         mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     (pow(voidfraction,1))
 );
#endif

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
 (
     IOobject
     (
        "phi",
        runTime.timeName(),
         mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(U*voidfraction) & mesh.Sf()
 );

Info<< "Generating interstitial face flux field phiByVoidfraction\n" << endl;
surfaceScalarField phiByVoidfraction
(
 IOobject
 (
    "phiByVoidfraction",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
 ),
 linearInterpolate(U) & mesh.Sf()
);

Info<< "Generating deltaY\n" << endl;
volScalarField deltaY
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0.*p
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

    surfaceScalarField alphacf("alphacf", fvc::interpolate(alpha));
    surfaceScalarField	alphaPhi("alphaPhi",alphacf*phi);

    const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();
forAll ( mesh.C(), celli)
{
    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
           pLocal[pointi] = pp[pLabels[pointi]];

    scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    deltaY[celli] = yDim;    // And similar for yDim and zDim
//      Info<<xDims<<endl;
}
singlePhaseTransportModel continiousPhaseTransport(U, phi);

   autoPtr<DPMIncompressibleTurbulenceModel<singlePhaseTransportModel> >
    turbulence
    (
        DPMIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
        (
            alpha,
            U,
            alphaPhi,
            phi,
            continiousPhaseTransport
        )
    );
