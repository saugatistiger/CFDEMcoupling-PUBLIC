
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
#include "transportModel.H"

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
#include <cmath>


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

        Info << "I ma here!!!!!!" << endl;
            #include "readTimeControls.H"
            #include "CourantNo.H"
       Info << "I am here!!!!!!!!!!!!!!!!!" << endl;
            #include "setDeltaT.H"
	////////////////////////////////////////////////////////////////////////////////////////
       	//volScalarField kField = turbulence->k();
	//volScalarField epsField = turbulence->epsilon();
	//volScalarField nutField = turbulence->nut();
	volScalarField magU = mag(U);

	//Info << "min(k): " << gMin(kField)
     	//<< " max(k): " << gMax(kField) << endl;

	//Info << "min(eps): " << gMin(epsField)
     	//	<< " max(eps): " << gMax(epsField) << endl;

	//Info << "max(nut): " << gMax(nutField) << endl;

	Info << "min(voidfraction): " << gMin(voidfraction) << endl;
	Info << "max(U): " << gMax(magU) << endl;
	////////////////////////////////////////////////////////////////////////////////////////
        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
	//bool hasEvolved = false;
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);
        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
        }

        Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
        Ksl.correctBoundaryConditions();
	voidfraction.correctBoundaryConditions();
	voidfraction = max(voidfraction, 0.1);
	voidfraction = min(voidfraction, 0.95);
        alpha = voidfraction;
	neeta = 1;
        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
        phi = voidfractionf*phiByVoidfraction;
	surfaceScalarField alphaf = fvc::interpolate(alpha);
        surfaceScalarField alphaPhi = phi*alphaf;
        gradAlpha = mag(fvc::grad(voidfraction));

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
	#include "deltaY_20260408.H"
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	forAll(alpha, i)
	{
		if (!std::isfinite(tau_vis[i].x()) ||
    		!std::isfinite(tau_vis[i].y()) ||
    		!std::isfinite(tau_vis[i].z()))
		{
    			Info << "Negative Negative!!!!!!" << endl;
		}
	}
        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
