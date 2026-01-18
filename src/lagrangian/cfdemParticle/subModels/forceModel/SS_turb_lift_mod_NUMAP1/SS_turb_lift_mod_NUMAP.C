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

#include "SS_turb_lift_mod_NUMAP.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>
//#include "mpi.h"
#include "fvMesh.H"
//#include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SS_turb_lift_mod_NUMAP, 0);

addToRunTimeSelectionTable
(
    forceModel,
    SS_turb_lift_mod_NUMAP,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SS_turb_lift_mod_NUMAP::SS_turb_lift_mod_NUMAP
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
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    gradAlpha_(sm.mesh().lookupObject<volScalarField>("gradAlpha")),
    useSecondOrderTerms_(false), sm_ (sm)
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

    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(false);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("lift_mod_NUMAPForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
    particleCloud_.probeM().vectorFields_.append("vorticity");  //other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");          //other are debug
    particleCloud_.probeM().scalarFields_.append("Rew");          //other are debug
    particleCloud_.probeM().scalarFields_.append("J_star");       //other are debug
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SS_turb_lift_mod_NUMAP::~SS_turb_lift_mod_NUMAP()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SS_turb_lift_mod_NUMAP::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();
    //treeSearch_(propsDict_.lookupOrDefault<Switch>("treeSearch", true));
    //searchEngine_(particleCloud_.mesh(),polyMesh::FACE_PLANES);
    scalar gradAlpha(0);
    vector position(0,0,0);
    vector position_up(0,0,0);
    vector position_dn(0,0,0);
    vector position1(0,0,0);
    scalar lift_mod_NUMAP(0);
    scalar L_turb_ratio(0);
    vector Us(0,0,0);
    vector Ur(0,0,0);
    vector U_up(0,0,0);
    vector U_dn(0,0,0);
    vector U_diff(0,0,0);
    vector Ur_mod(0,0,0);
    scalar magU_diff(0);
    scalar magUr(0);
    scalar magVorticity(0);
    scalar ds(0);
    scalar U_star(0);
    scalar U_star_new(0);
    scalar dParcel(0);
    scalar dragCoefficient(0);
    vector drag(0,0,0);
    scalar drag_y(0);
    scalar nuf(0);
    scalar Ks_st(0);
    scalar rho(0);
    scalar voidfraction(1);
    scalar Rep(0);
    scalar Rew(0);
    scalar Cl(0);
    scalar Cl1(0);
    scalar Cl2_d(0);
    scalar Cd(0);
    scalar Cl_star(0);
    scalar J_star(0);
    scalar Omega_eq(0);
    scalar alphaStar(0);
    scalar epsilon(0);
    scalar omega_star(0);
    scalar H_ = readScalar(propsDict_.lookup("H"));
    vector vorticity(0,0,0);
    volVectorField vorticity_ = fvc::curl(U_);  
    volTensorField Ugrad_all = fvc::grad(U_); 
    vector Ur_u(0,0,0);
    vector Ur_l(0,0,0);
    vector force(0,0,0);
    vector gradU(0,0,0);
    vector gradU_int(0,0,0);
    vector gradU_mod(0,0,0);
    vector gradU_mod11(0,0,0);
    vector ab(0,0,0);
    scalar tau(0);

    #include "resetVorticityInterpolator.H"
    #include "resetUInterpolator.H"

    #include "setupProbeModel.H"

    const faceList & ff = sm_.mesh().faces();
    const pointField & pp = sm_.mesh().points(); 

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            force          = vector::zero;
            int cellI    = particleCloud_.cellIDs()[index][0];
            position       = particleCloud_.position(index);

            int foundIndex = -1;
            for (int i = 0; i < particleCloud_.numberOfParticles(); ++i) 
            {
            if (particleCloud_.cellIDs()[i][0] == cellI) 
                {
                   foundIndex = i;
                   break; 
                }
            }
            vector position1 = particleCloud_.position(index);
            if (cellI > -1) // particle Found
            {
                nuf = nufField[cellI];
                rho = rhoField[cellI];
                Us  = particleCloud_.velocity(index);

                ds  = 2*particleCloud_.radius(index);
                if( forceSubM(0).interpolation() )
                {
                    Ur             = UInterpolator_().interpolate(position,cellI) - Us;
                    vorticity      = vorticityInterpolator_().interpolate(position,cellI);
                    Ur_u           = UInterpolator_().interpolate(position + vector (ds/2,ds/2,ds/2),cellI) - Us;
                    Ur_l           = UInterpolator_().interpolate(position - vector (ds/2,ds/2,ds/2),cellI) - Us;
                    gradU          = (Ur_u-Ur_l)/ds;
                    gradU_mod	   = vector (fabs(gradU[1]),fabs(gradU[0]),fabs(gradU[2]));
                }
                else
                {
                    Ur =  UInterpolator_().interpolate(position,cellI) - Us;
                    vorticity = vorticity_[cellI];
  		    gradU_mod	   = vector (sqrt(fabs(Ugrad_all[cellI][1])),sqrt(fabs(Ugrad_all[cellI][3])),sqrt(fabs(Ugrad_all[cellI][8])));
                }
                U_diff = U_up - U_dn;
                magU_diff = mag(magU_diff);
                gradAlpha = gradAlpha_[cellI];
                magUr           = mag(Ur);
                if (magUr > 0)
                {
                     // calc fluid drag Coeff
                    Rep = ds*voidfraction*magUr/(nuf+SMALL);
                    Cd = sqr(0.63 + 4.8/sqrt(Rep));
 		    scalar Xi = 3.7 - 0.65 * exp(-sqr(1.5-log10(Rep))/2);
                    dragCoefficient = 0.125*Cd*rho*M_PI*ds*ds*pow(voidfraction,(2-Xi))*magUr;
                    if (modelType_=="B")
                        dragCoefficient /= voidfraction;

                    drag = dragCoefficient*Ur;
                    
                    scalar gradAlpha_ch = gradAlpha/(*std::max_element(gradAlpha_.begin(), gradAlpha_.end()));
                    if (gradAlpha_ch > 0.5)
                        {
                        tau = (nuf*fabs(Ugrad_all[cellI][3]));
			//Pout << "Entered gradAlpha_ch" << endl;
                        }
                    else
                        {
                        tau = (0);
                        }
                    const cell & cc = sm_.mesh().cells()[cellI];
                    labelList pLabels(cc.labels(ff));
                    pointField pLocal(pLabels.size(), vector::zero);

                    forAll (pLabels, pointi)
                        pLocal[pointi] = pp[pLabels[pointi]];

                    scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
                    U_star      = sqrt(tau/rho);//sm_.mesh().points()[cellI_up]
                    scalar Yplus   = 2*yDim*rho*U_star/nuf;
                    scalar Kspl = 100;//rho*2.5*ds*U_star/nuf;
                    drag = dragCoefficient*Ur;
                    if (gradAlpha_ch > 0.5)
                    {

				////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////// NUMAP  ////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////
           	 		Foam::vector position_up (position[0], position[1]+ds, position[2]);
            			Foam::vector position_dn (position[0], position[1]-ds, position[2]);
            			int cellI_up   = searchEngine_.findCell(position_up,-1,treeSearch_);
            			int cellI_dn   = searchEngine_.findCell(position_dn,-1,treeSearch_);
	    			U_up = U_[cellI_up];
				if (cellI_dn > -1)
					{
            					U_dn = U_[cellI_dn];
					}
				else
					{
						U_dn = vector(0,0,0);
					}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////
                        U_star      = (0.001*sqrt(tau/rho));
		                        if (Yplus >= 5)
                        {
			
                            for (int i = 0; i < 40; ++i) 
                            {
				//Pout << "U_star" << U_star << endl;
	
                                scalar term1 = (1 / 0.41) * std::log((9.81 * yDim * U_star) / nuf);
				//Pout << "term1 = "<< term1 << endl;
                                scalar term2 = (abs(U_diff[0]) / U_star);
				//Pout << "term2 = "<< term2 << endl;
                                scalar term3 = 0.0;
                                scalar term4 = 0.0;
				//Pout << "22222222222222222" << endl;
                                if (Kspl>=90)
                                {
                                term3 = (1 / 0.41) * std::log(1+((0.5 * ds*2.5 * U_star) / nuf));
				//Pout << "term3 = "<< term3 << endl;
                                term4 = ((0.5*ds*2.5*rho)/(0.5*2.5*ds*rho*U_star+rho*nuf));
				//Pout << "term4 = "<< term4 << endl;
                                }
                                else if (Kspl>=2.25 && Kspl<90)
                                {
                                term3 = (1 / 0.41) *((((rho*ds*2.5*U_star/nuf)-2.25)/87.75)+(0.5*rho*2.5*ds*U_star/nuf))*std::sin(0.4258*(std::log(rho*2.5*ds*U_star/nuf)-0.811));

                                scalar A = (1/0.41)*((((0.5*2.5*ds*rho)/nuf)) + (4*2.5*ds*rho/(351*nuf)))*std::sin(0.4258*(std::log(rho*2.5*ds*U_star/nuf)-0.811));
                                scalar B = (2129/(5000*0.41*U_star))*((4/351)*((rho*2.5*ds*U_star/nuf)-2.25)+(0.5*2.5*ds*rho/nuf))*std::cos(0.4258*(std::log(rho*2.5*ds*U_star/nuf)-0.811));
                                term4 = A + B;
                                }
                                else if (Kspl<2.25)
                                {
                                term3 = 0;
                                term4 = 0;
                                }
                                scalar numerator = term1 - term2 - term3;
                                scalar denominator = (1 / (0.41 * U_star)) + (abs(U_diff[0]) / (U_star * U_star)) - term4;
				
                                U_star_new = U_star - (numerator / denominator); 
                                U_star = (U_star_new);
   				//Pout << "numerator = " << numerator << endl;
				//Pout << "denominator = " <<denominator << endl;
 				
                            }
                        drag_y   = (2.0/8.0)*U_star*U_star*3.14*rho*ds*ds;
                        }
                        else if (Yplus < 5)
                        {
                            drag_y   = (2.0/8.0)*tau*3.14*ds*ds;
                        }
                    }
		   
                    drag[0]  = drag_y/voidfraction  + drag[0];
                    Ks_st = U_star*ds/(nuf);
                    //if (gradAlpha_ch > 0.5) {
                    if (U_star != 0)
                    {
                        if (Ks_st<=60)
                        {
                            L_turb_ratio = 3 * std::exp(-std::pow(log10(Ks_st) - 1.32, 2) / (2 * std::pow(0.41, 2)));
                        }
                        else
                        {
                            L_turb_ratio = 1.9;
                        }
                    }
                    else
                    {
                        L_turb_ratio = 0;
                    }
		    
                    force[0] = drag[1]*L_turb_ratio + drag[0];
		            force[1] = drag[0]*L_turb_ratio + drag[1];
                    force[2] = drag[2];

                    forceSubM(0).scaleForce(force,dParcel,index);

                }

                //**********************************        
                //SAMPLING AND VERBOSE OUTOUT
                if( forceSubM(0).verbose() )
                {   
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "dprim = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
    		    Pout << "lift_mod_NUMAP = " << lift_mod_NUMAP << endl;
                Pout << "H = " << H_ << endl;
                }
		
                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, force);           //first entry must the be the force
                    vValues.setSize(vValues.size()+1, Ur);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
                // END OF SAMPLING AND VERBOSE OUTOUT
                //**********************************        

            }
            // write particle based data to global array
	    
            forceSubM(0).partToArray(index,force,vector::zero);
        //}
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
