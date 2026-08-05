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

     This function is based on the derivation in R. Kurose, S. Komori, 
     Drag and lift_mod_20251111 forces on a rotating sphere in a linear shear flow, 
     Journal of Fluid Mechanics. 384 (1999) 183-206.

     The data for this functions is based on J.B. Mclaughlin, 
     Inertial migration of a small sphere in linear shear flows, 
     Journal of Fluid Mechanics. 224 (1991) 261-274.

     The second order terms are based on: 
     Mei lift_mod_20251111 force following Loth and Dorgan 2009,
     Environ Fluid Mech (2009) 9:187–206, DOI 10.1007/s10652-009-9123-x
     and can be added to the lift_mod_20251111 coefficient if desired
     (contribution from RQ)

    - including interpolation of the velocity to the particle position (optional)
    - including output to file for testing/data analysis (optional)


Class
    SS_turb_no_lift_mod_20260603

SourceFiles
    SS_turb_no_lift_mod_20260603.C

\*---------------------------------------------------------------------------*/

#ifndef SS_turb_no_lift_mod_20260603_H
#define SS_turb_no_lift_mod_20260603_H

#include "forceModel.H"
#include "IOList.H"
#include "interpolationCellPoint.H"
# include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SS_turb_no_lift_mod_20260603 Declaration
\*---------------------------------------------------------------------------*/

class SS_turb_no_lift_mod_20260603
:
    public forceModel
{
protected:

    dictionary propsDict_;

    Switch treeSearch_;

    mutable meshSearch searchEngine_;
    
    //const Foam::fvMesh& mesh_;



private:
    //dictionary propsDict_;

   cfdemCloud& sm_;

    word velFieldName_;

    word velFieldName1_;

    word H_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    //const volScalarField& k_;

    const volVectorField& LiftDir_;

    const volVectorField& tau_vis_;

    bool useSecondOrderTerms_;

    word UsFieldName_;   

    const volVectorField& UsField_; 

public:

    //- Runtime type information
    TypeName("SS_turb_no_lift_mod_20260603");


    // Constructors

        //- Construct from components
        SS_turb_no_lift_mod_20260603
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word name = word("")
        );

    // Destructor

        ~SS_turb_no_lift_mod_20260603();


    // Member Functions
        void setForce() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
