
/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

#include "OFversion.H"
#include "DPMIncompressibleTurbulenceModel.H"
#include "pisoControl.H"
    #include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "cfdemCloud.H"

#if defined(anisotropicRotation)
    #include "cfdemCloudRotation.H"
#endif
#if defined(superquadrics_flag)
    #include "cfdemCloudRotationSuperquadric.H"
#endif
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
        pisoControl piso(mesh);
        #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAccleration_new.H"
    #include "checkImCoupleM.H"
    #if defined(anisotropicRotation)
        cfdemCloudRotation particleCloud(mesh);
    #elif defined(superquadrics_flag)
        cfdemCloudRotationSuperquadric particleCloud(mesh);
    #else
        cfdemCloud particleCloud(mesh);
    #endif
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"

        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        Info<< "I am Here" << endl;
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);
        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
        }
    
        Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
        Ksl.correctBoundaryConditions();
	voidfraction.correctBoundaryConditions();
        alpha = (pow(voidfraction,1));
        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
        phi = voidfractionf*phiByVoidfraction;
	surfaceScalarField alphaf = fvc::interpolate(alpha);
        surfaceScalarField alphaPhi = phi*alphaf;
        gradAlpha = mag(fvc::grad(voidfraction));
//	forAll(owner, facei)
// 	{
// 	delta[facei] = C[neighbour[facei]] - C[owner[facei]];
// 	}
	#include "deltaY.H"
///////////////////////////////////////////////////////////////////////////////        //Force Checks
        #include "forceCheckIm.H"

        //#include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            {
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
                  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
//                + turbulence->divDevReff(U)
                  + particleCloud.divVoidfractionTau(U, voidfraction)
                  ==
                  - fvm::Sp(Ksl/rho,U)
                  + fvOptions(U)
                );

                UEqn.relax();
                fvOptions.constrain(UEqn);

                    if (piso.momentumPredictor())
                {
                    if (modelType=="B" || modelType=="Bfull")
                        solve(UEqn == - fvc::grad(p) + Ksl/rho*Us);
                    else
                        solve(UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us);

                    fvOptions.correct(U);
                }

                // --- PISO loop
                    while (piso.correct())
                {
                    volScalarField rUA = 1.0/UEqn.A();

                    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
                    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));

                    U = rUA*UEqn.H();
                        phi = ( fvc::interpolate(U) & mesh.Sf() )
                            + rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction);

                    surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
                    phi += rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                    if (modelType=="A")
                        rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                    // Update the fixedFluxPressure BCs to ensure flux consistency
                    #include "fixedFluxPressureHandling.H"
                    
                    // Non-orthogonal pressure corrector loop
  		    while (piso.correctNonOrthogonal())
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rUAvoidfraction, p) == fvc::div(voidfractionf*phi) + particleCloud.ddtVoidfraction()
                        );
                        pEqn.setReference(pRefCell, pRefValue);

			pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                            if (piso.finalNonOrthogonalIter())
                            {
                                phiByVoidfraction = phi - pEqn.flux()/voidfractionf;
                            }
                     

                    } // end non-orthogonal corrector loop

                    phi = voidfractionf*phiByVoidfraction;
                    //alphaPhi = voidfractionf;
                    #include "continuityErrorPhiPU.H"
              
		    	
                    if (modelType=="B" || modelType=="Bfull")
                        U -= rUA*fvc::grad(p) - Ksl/rho*Us*rUA;
                    else
                        U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA;

                    U.correctBoundaryConditions();
                    fvOptions.correct(U);

                } // end piso loop
            }

            continiousPhaseTransport.correct();
            turbulence->correct();
        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
