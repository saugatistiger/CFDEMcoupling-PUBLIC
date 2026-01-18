/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "SS_turb_lift_mod_20251111.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>
//#include "mpi.h"
#include "fvMesh.H"
//#include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SS_turb_lift_mod_20251111, 0);

addToRunTimeSelectionTable
(
    forceModel,
    SS_turb_lift_mod_20251111,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SS_turb_lift_mod_20251111::SS_turb_lift_mod_20251111
(
    const dictionary& dict,
    cfdemCloud& sm,
    word name
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    treeSearch_(propsDict_.lookupOrDefault<Switch>("treeSearch", true)),
    searchEngine_(particleCloud_.mesh(),polyMesh::FACE_PLANES),
    velFieldName_(propsDict_.lookup("velFieldName")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    tau_vis_(sm.mesh().lookupObject<volVectorField>("tau_vis")),
    LiftDir_(sm.mesh().lookupObject<volVectorField>("LiftDir")),
    k_(sm.mesh().lookupObject<volScalarField>("k")),
    useSecondOrderTerms_(false), sm_ (sm),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_))
{

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch
    //propsDict_(dict.subDict(name == "" ? typeName + "Props" : name + "Props"));
    //set default switches (hard-coded default = false)
    forceSubM(0).setSwitches(0,true);  // enable treatExplicit, otherwise this force would be implicit in slip vel! - IMPORTANT!
    forceSubM(0).setSwitches(1,true); 

    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(false);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
    particleCloud_.probeM().vectorFields_.append("vorticity");  //other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");          //other are debug
    particleCloud_.probeM().scalarFields_.append("Rew");          //other are debug
    particleCloud_.probeM().scalarFields_.append("J_star");       //other are debug
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SS_turb_lift_mod_20251111::~SS_turb_lift_mod_20251111()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SS_turb_lift_mod_20251111::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();
    //treeSearch_(propsDict_.lookupOrDefault<Switch>("treeSearch", true));
    //searchEngine_(particleCloud_.mesh(),polyMesh::FACE_PLANES);
    vector tau_vis(0,0,0);
    vector LiftDir(0,0,0);
    vector position(0,0,0);
    scalar lift_mod_20251111(0);
    vector dragExplicit(0,0,0);
    vector Ufluid(0,0,0);
    scalar L_turb_ratio(0);
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar magUr(0);
    scalar ds(0);
    scalar U_star(0);
    scalar dParcel(0);
    scalar dragCoefficient(0);
    vector drag(0,0,0);
    scalar nuf(0);
    scalar rho(0);
    scalar voidfraction(1);
    scalar k(0);
    scalar Rep(0);
    scalar Rew(0);
    scalar Cd(0);
    vector vorticity(0,0,0);
    volVectorField vorticity_ = fvc::curl(U_);  
    vector force(0,0,0);

    #include "resetVorticityInterpolator.H"
    #include "resetUInterpolator.H"
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            force          = vector::zero;
            int cellI      = particleCloud_.cellIDs()[index][0];
            position       = particleCloud_.position(index);
            if (cellI > -1) // particle Found
            {
                nuf = nufField[cellI];
                rho = rhoField[cellI];
                Us  = particleCloud_.velocity(index);
                ds  = 2*particleCloud_.radius(index);
                position       = particleCloud_.position(index);
                if( forceSubM(0).interpolation() )
                {
                    voidfraction   = voidfractionInterpolator_().interpolate(position,cellI);
                    Ufluid         = UInterpolator_().interpolate(position,cellI);
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.30) voidfraction = 0.30;
                }
                else
                {
                    Ufluid = U_[cellI];
                    voidfraction = voidfraction_[cellI];
                }
                k = k_[cellI];
		        tau_vis	= tau_vis_[cellI];
                LiftDir = LiftDir_[cellI];
                Ur = Ufluid-Us;
                //Info << "Ur_before = " << Ur << endl;
                
                int dir = 1;
                /*
                if (mag(tau_vis)>0)   
                {
                    for (int ii=0;ii<3;++ii)
                    {
                        if ((Ur[ii]!=0) && (mag(Ur)!=0))
                        {    
                            if (Ur[ii]<0)
                            {
                                dir = -1;
                            }
                            Ur[ii] = dir*Foam::sqrt(Ur[ii]*Ur[ii] + (abs(Ur[ii])/mag(Ur))*0.6*k);
                        }
                    }
                }
                */
                //Info << "Ur_after = " << Ur << endl;
                magUr   = mag(Ur);
                if (magUr > 0)
                {
                    Rep = ds*voidfraction*magUr/(nuf+SMALL);
                    Cd = sqr(0.63 + 4.8/sqrt(Rep));
 		            scalar Xi = 3.7 - 0.65 * exp(-sqr(1.5-log10(Rep))/2);
                    dragCoefficient = 0.125*Cd*rho*M_PI*ds*ds*pow(voidfraction,(2-Xi))*magUr;
                    dragCoefficient /= voidfraction;
                    drag = dragCoefficient*Ur;
		            tau_vis   = tau_vis_[cellI];
		            vector drag_vis = (tau_vis*ds*ds*M_PI*0.25)/voidfraction;
                    drag  = drag_vis  + drag;
                    vector lift(0,0,0);
                    vector L_turb_ratio(0,0,0);
                    vector Ks_st(0,0,0);
                    vector U_star(0,0,0);
                    vector drag_vis_lift(0,0,0);
                    U_star[0] = Foam::sqrt((mag(tau_vis[1])+mag(tau_vis[2]))/rho);
                    Ks_st[0] = U_star[0]*ds/(nuf);

                    if (mag(drag_vis)>0)
                    {
                        drag_vis_lift = drag;
                    }
                    
                    if (U_star[0] != 0)
                    {
                        L_turb_ratio[0] = 2.5*std::exp(-std::pow(log10(Ks_st[0])-1.25,2)/0.18) + 1.5*(1-std::exp(-0.1*Ks_st[0]*Ks_st[0]));   
                    }
                    lift[0] = (((Foam::sqrt(std::pow(drag_vis_lift[1],2)+std::pow(drag_vis_lift[2],2)))  * L_turb_ratio[0]))*LiftDir[0];

                    U_star[1] = Foam::sqrt((mag(tau_vis[0])+mag(tau_vis[2]))/rho);
                    Ks_st[1] = U_star[1]*ds/(nuf);
                    if (U_star[1] != 0)
                    {
                        L_turb_ratio[1] = 2.5*std::exp(-std::pow(log10(Ks_st[1])-1.25,2)/0.18) + 1.5*(1-std::exp(-0.1*Ks_st[1]*Ks_st[1]));   
                    }
                    lift[1] = (((Foam::sqrt(std::pow(drag_vis_lift[0],2)+std::pow(drag_vis_lift[2],2)))* L_turb_ratio[1]))*LiftDir[1];

                    U_star[2] = Foam::sqrt((mag(tau_vis[0])+mag(tau_vis[1]))/rho);
                    Ks_st[2] = U_star[2]*ds/(nuf);
                    if (U_star[2] != 0)
                    {
                        L_turb_ratio[2] = 2.5*std::exp(-std::pow(log10(Ks_st[2])-1.25,2)/0.18) + 1.5*(1-std::exp(-0.1*Ks_st[2]*Ks_st[2]));   
                    }
                    lift[2] = (((Foam::sqrt(std::pow(drag_vis_lift[0],2)+std::pow(drag_vis_lift[1],2)))* L_turb_ratio[2]))*LiftDir[2];
                    force = lift + drag_vis;
                    //force[0] = L_turb_ratio[0]*LiftDir[0]*0.0001;//lift + drag;
                    //force[1] = L_turb_ratio[1]*LiftDir[1]*0.0001;
                    //force[2] = L_turb_ratio[2]*LiftDir[2]*0.0001;
                    /*
                    Info << "           "<< endl;
                    Info << "LiftDir = "<<LiftDir << endl;
                    Info << "U_star = "<<U_star << endl;
                    Info << "Ks_st = "<<Ks_st << endl;
                    Info << "L_turb_ratio = "<<L_turb_ratio << endl;
                    Info << "Drag_vis = "<<drag_vis << endl;
                    Info << "Lift = "<<lift << endl;
                    Info << "Drag_total = "<<drag << endl;
                    */
                    //forceSubM(0).scaleCoeff(dragCoefficient,dParcel,index);
                    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                        forceSubM(iFSub).explicitCorr( drag, 
                                                       dragExplicit,
                                                       dragCoefficient,
                                                       Ufluid, U_[cellI], Us, UsField_[cellI],
                                                       forceSubM(iFSub).verbose()
                                                     );
                }
                if( forceSubM(0).verbose() )
                {   
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "dprim = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
    		    Pout << "lift_mod_20251111 = " << lift_mod_20251111 << endl;
                }
		
                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.setSize(vValues.size()+1, force);
                    vValues.setSize(vValues.size()+1, Ur);
                    sValues.setSize(sValues.size()+1, Rep); 
                    sValues.setSize(sValues.size()+1, Cd);
                    sValues.setSize(sValues.size()+1, voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
            // write particle based data to global array    
            forceSubM(0).partToArray(index,force,dragExplicit,Ufluid,dragCoefficient);
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

