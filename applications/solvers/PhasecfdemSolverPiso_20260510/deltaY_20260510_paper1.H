
neeta = dimensionedScalar("zero", dimless, 0.0);
const volScalarField& nu = mesh.lookupObject<volScalarField>("nu");
forAll ( alpha, i)
{
	if (alpha[i] < 0.76)
	{
    		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			if (alpha[adj[j]]<0.76)
			{
				neeta[i] = 1;
				break;
			}
		}
		scalar n = 0;
		for (int j=0; j < adj.size(); ++j)
		{
			if (alpha[adj[j]]<0.76)// && mag(Us[adj[j]])<u_star.value())
			{
				n = n+1;
			}
		}
		if (n == adj.size())
		{
			neeta[i] = 0;
		}
		
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

//forAll ( alpha, i)
//{
//	if (neeta[i] > 0)
//	{
//		scalar nnn = 0;
//		labelList adj = mesh.cellCells()[i];
//		for (int j = 0; j < adj.size(); ++j)
//		{
//			if (neeta[adj[j]] > 0)
//			{
//				nnn = nnn + 1;
//			}
//		}
//		if (nnn==0)
//		{
//			neeta[i] = 0;
//		}
//	}
//}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
forAll ( alpha, i)
{
	if (((neeta[i] > 0) && (mag(Us[i][0])>4*u_star.value()) && mesh.C()[i][1]>4*D50.value()))
	{
		neeta[i] = 0;
	}
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// =========================================================================
//  Step 5 — Remove isolated neeta cells (no neeta neighbours)
// =========================================================================

forAll(alpha, i)
{
    if (neeta[i] > 0)
    {
        const labelList& adj = mesh.cellCells()[i];
        scalar nnn = 0;
        for (int j = 0; j < adj.size(); ++j)
        {
            if (neeta[adj[j]] > 0) nnn++;
        }
        if (nnn == 0)
        {
            neeta[i] = 0;
        }
    }
}

// =========================================================================
//  Step 3 — Spread neeta horizontally in x (4 passes)
// =========================================================================
/*
for (int pass = 0; pass < 4; ++pass)
{
    forAll(alpha, i)
    {
        if (neeta[i] > 0)
        {
            const labelList& adj = mesh.cellCells()[i];
            for (int j = 0; j < adj.size(); ++j)
            {
                vector D = mesh.C()[adj[j]] - mesh.C()[i];
                if ((mag(D[0]) > 0.001) && (D[1] < 0.001) && (mag(D[2]) < 0.001))
                {
                    neeta[adj[j]] = 1;
                }
            }
        }
    }
}

// =========================================================================
//  Step 4 — Spread neeta vertically in z (4 passes)
// =========================================================================
for (int pass = 0; pass < 4; ++pass)
{
    forAll(alpha, i)
    {
        if (neeta[i] > 0)
        {
            const labelList& adj = mesh.cellCells()[i];
            for (int j = 0; j < adj.size(); ++j)
            {
                vector D = mesh.C()[adj[j]] - mesh.C()[i];
                if ((mag(D[0]) < 0.001) && (D[1] < 0.001) && (mag(D[2]) > 0.001))
                {
                    neeta[adj[j]] = 1;
                }
            }
        }
    }
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// =========================================================================
//  Step 6 — Ensure exactly one neeta=1 cell per x-z column (lowest y)
//            and gap-fill empty columns using BFS + alpha gradient prediction
// =========================================================================
{
    // Compute gradAlpha globally
    volVectorField gradAlpha = fvc::grad(alpha);

    // Compute x-z cell spacing directly from mesh neighbour distances
    scalar minDx = GREAT;
    scalar minDz = GREAT;

    forAll(alpha, i)
    {
        const labelList& adj = mesh.cellCells()[i];
        for (int j = 0; j < adj.size(); ++j)
        {
            vector D = mesh.C()[adj[j]] - mesh.C()[i];

            if ((mag(D[0]) > SMALL) && (mag(D[1]) < mag(D[0])) && (mag(D[2]) < mag(D[0])))
            {
                if (mag(D[0]) < minDx) minDx = mag(D[0]);
            }

            if ((mag(D[2]) > SMALL) && (mag(D[1]) < mag(D[2])) && (mag(D[0]) < mag(D[2])))
            {
                if (mag(D[2]) < minDz) minDz = mag(D[2]);
            }
        }
    }

    const scalar xzTol = min(minDx, minDz) * 0.5;

    auto bucketCoord = [&](scalar val) -> int {
        return static_cast<int>(std::round(val / xzTol));
    };

    // columnBed      : (ix,iz) -> (lowest_y, cellIndex) among neeta=1 cells
    // columnAllCells : (ix,iz) -> all cells (y, cellIndex) in that column
    std::map<std::pair<int,int>, std::pair<scalar,label>>              columnBed;
    std::map<std::pair<int,int>, std::vector<std::pair<scalar,label>>> columnAllCells;

    forAll(alpha, i)
    {
        const vector& C = mesh.C()[i];
        auto key = std::make_pair(bucketCoord(C[0]), bucketCoord(C[2]));

        columnAllCells[key].push_back({C[1], i});

        if (neeta[i] > 0)
        {
            auto it = columnBed.find(key);
            if (it == columnBed.end() || C[1] < it->second.first)
            {
                columnBed[key] = {C[1], i};
            }
        }
    }

    // Wipe neeta and re-mark only the single lowest cell per column
    neeta = dimensionedScalar("zero", dimless, 0.0);
    for (const auto& entry : columnBed)
    {
        neeta[entry.second.second] = 1;
    }

    // Breadth-first gap-fill seeded from all already-filled columns
    std::queue<std::pair<int,int>> fillQueue;
    for (const auto& entry : columnBed)
    {
        fillQueue.push(entry.first);
    }

    while (!fillQueue.empty())
    {
        auto key = fillQueue.front();
        fillQueue.pop();

        label  srcCell = columnBed[key].second;
        scalar srcY    = columnBed[key].first;
        vector gA      = gradAlpha[srcCell];

        for (int dx = -1; dx <= 1; ++dx)
        for (int dz = -1; dz <= 1; ++dz)
        {
            if (dx == 0 && dz == 0) continue;

            auto nKey = std::make_pair(key.first + dx, key.second + dz);

            if (columnBed.count(nKey))       continue;
            if (!columnAllCells.count(nKey)) continue;

            // Predict bed height using alpha gradient iso-surface condition
            scalar dispX = dx * xzTol * 2;
            scalar dispZ = dz * xzTol * 2;
            scalar predictedY = srcY;
            if (mag(gA.y()) > SMALL)
            {
                predictedY = srcY - (gA.x()*dispX + gA.z()*dispZ) / gA.y();
            }

            label  bestCell  = -1;
            scalar bestDist  = GREAT;

            // Tier 1: alpha < 0.76 closest to predictedY
            for (const auto& cellEntry : columnAllCells[nKey])
            {
                if (alpha[cellEntry.second] < 0.76)
                {
                    scalar dist = mag(cellEntry.first - predictedY);
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestCell = cellEntry.second;
                    }
                }
            }

            // Tier 2: alpha < 0.9 closest to predictedY
            if (bestCell < 0)
            {
                bestDist = GREAT;
                for (const auto& cellEntry : columnAllCells[nKey])
                {
                    if (alpha[cellEntry.second] < 0.9)
                    {
                        scalar dist = mag(cellEntry.first - predictedY);
                        if (dist < bestDist)
                        {
                            bestDist = dist;
                            bestCell = cellEntry.second;
                        }
                    }
                }
            }

            // Tier 3: no alpha filter — just closest to predictedY
            if (bestCell < 0)
            {
                bestDist = GREAT;
                for (const auto& cellEntry : columnAllCells[nKey])
                {
                    scalar dist = mag(cellEntry.first - predictedY);
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestCell = cellEntry.second;
                    }
                }
            }

            if (bestCell >= 0)
            {
                neeta[bestCell] = 1;
                columnBed[nKey] = {mesh.C()[bestCell][1], bestCell};
                fillQueue.push(nKey);
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
scalarList 	innerList(6,0.0);
scalarListList 	neeta_nei(alpha.size(),innerList);
scalarListList 	neeta_nei_alpha(alpha.size(),innerList);
scalarListList 	neeta_dir(alpha.size(),innerList); // neighbour of the neeta
vector DD (0,0,0);
forAll ( alpha, i)
{
	int n = 0;
	if ( neeta[i] == 1 )
	{
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j<adj.size();++j)
		{
			if (alpha[adj[j]] > 0.6)
			{
				//Info << "n=" << n << endl;
				//Info << "alpha[adj[j]]" << alpha[adj[j]] << endl;
				neeta_nei[i][n]	=	adj[j];
				//
				DD = mesh.C()[i] - mesh.C()[neeta_nei[i][n]];
				//Info << "DD =" << DD << endl;
				if ((mag(DD[0]) > mag(DD[1])) && (mag(DD[0]) > mag(DD[2])))
				{
					if (DD[0]>0)
					{
						neeta_dir[i][n] = 1;
					}
					else
					{
						neeta_dir[i][n] = -1;
					}
				}
				else if ((mag(DD[1]) > mag(DD[0])) && (mag(DD[1]) > mag(DD[2])))
				{
					if (DD[1]>0)
					{
						neeta_dir[i][n] = 2;
					}
					else
					{
						neeta_dir[i][n] = -2;
					}
				}
				else if ((mag(DD[2]) > mag(DD[0])) && (mag(DD[2]) > mag(DD[1])))
				{
					if (DD[2]>0)
					{
						neeta_dir[i][n] = 3;
					}
					else
					{
						neeta_dir[i][n] = -3;
					}
				}
				neeta_nei_alpha[i][n] = (alpha[adj[j]]);
				//Info << "neeta_nei[i][n]="<<alpha[neeta_nei[i][n]] << endl;
				//Info << "adj[j]="<<alpha[adj[j]] << endl;
				//Info << "neeta_nei_alpha[i][j] =" <<neeta_nei_alpha[i] << endl;
				//
				n = n +1;
			}
		}
		//Info << "neeta_nei[i]="<<neeta_nei[i] << endl;
		//Info << "neeta_dir[i]="<<neeta_dir[i] << endl;
		//Info << "neeta_nei_alpha[i][j] =" <<neeta_nei_alpha[i] << endl;
		for (int jj =0; jj < 3; ++jj)
		{
			int nn=0;
			scalarList  Di_i(2,0.0);
			scalarList  Di(2,0.0);
			for (int ii =0; ii < 6; ++ii)
			{
				if (mag(neeta_dir[i][ii])==jj+1)
				{
					Di_i[nn] = neeta_nei[i][ii];
					Di[nn] = neeta_dir[i][ii];
					nn = nn+1;
				}
			}
			//Info << "I am here!!!!!!!!!!!" << endl;
			//Info << Di << endl;
			if (Di[0]==0)
			{
				LiftDir[i][jj] = 0;
			}
			else if (Di[1]!=0)
			{
				if (alpha[Di_i[0]]> alpha[Di_i[1]])
				{
					LiftDir[i][jj] = - Di[0]/mag(Di[0]);
				}
				else
				{
					LiftDir[i][jj] = - Di[1]/mag(Di[1]);
				}
			}
			else
			{
				LiftDir[i][jj] = - Di[0]/mag(Di[0]);
			}
			//Info << "I am  not here!!!!!!!!!!!" << endl;
		}
	}
}

forAll (alpha, i)
{
	vector U_star_dir(0,0,0);
	if ( neeta[i] == 1 )
	{
	int a = 0;
	for (int k = 0; k < 6 ; ++k)
	{
		if (neeta_nei[i][k]==0)
		{
			break;
		}
		else
		{
			a = a+1;	
		}
	}
		vector D (0,0,0);
		vector U_rel (0,0,0);
		vector U_plane (0,0,0);
		vector alpha_angle(0,0,0);
		vector alpha_angle_res(0,0,0);
		vector beta_angle(0,0,0);
		vector beta_angle_res(0,0,0);
		vector U_star_plane(0,0,0);
		vector Yplus(0,0,0);
		scalar AA(0);
		scalar BB(0);
		for (int j=0 ; j < a ; ++j)
		{
			D = mesh.C()[neeta_nei[i][j]] - mesh.C()[i];
			for (int ii=0;ii<3;++ii)
            {
                if ((U[neeta_nei[i][j]][ii]!=0) && (mag(U[neeta_nei[i][j]])!=0))
                {
                	//AA = Foam::sqrt(U[neeta_nei[i][j]][ii]*U[neeta_nei[i][j]][ii] + (mag(U[neeta_nei[i][j]][ii])/mag(U[neeta_nei[i][j]]))*0.6*k[neeta_nei[i][j]]);
			AA = Foam::sqrt(U[neeta_nei[i][j]][ii]*U[neeta_nei[i][j]][ii] + (mag(U[neeta_nei[i][j]][ii])/mag(U[neeta_nei[i][j]]))*0.6*0);//k[neeta_nei[i][j]]);
				}
				if ((U[i][ii]!=0) && (mag(U[i])!=0))
				{
					//BB = Foam::sqrt(U[i][ii]*U[i][ii] + (mag(U[i][ii])/mag(U[i]))*0.6*k[i]);
					BB = Foam::sqrt(U[i][ii]*U[i][ii] + (mag(U[i][ii])/mag(U[i]))*0.6*0);//k[i]);
				}
				U_rel[ii] = AA - BB;
			}
			
			#define SAFE_ACOS(x) Foam::acos(Foam::max(Foam::min(scalar(x), scalar(1.0)), scalar(-1.0)))
			
			U_plane[0] = Foam::sqrt(U_rel[1]*U_rel[1] + U_rel[2]*U_rel[2]); //yz-plane
			alpha_angle[0] = SAFE_ACOS(mag(U_rel[1])/(mag(U_plane[0])+SMALL));
			beta_angle[0]  = SAFE_ACOS(mag(U_rel[2])/(mag(U_plane[0])+SMALL));
			alpha_angle_res[0] = Foam::cos(alpha_angle[0]); //
			beta_angle_res[0] = Foam::cos(beta_angle[0]);
			
			U_plane[1] = Foam::sqrt(U_rel[0]*U_rel[0] + U_rel[2]*U_rel[2]);
			alpha_angle[1] = SAFE_ACOS(mag(U_rel[0])/(mag(U_plane[1])+SMALL));
			beta_angle[1]  = SAFE_ACOS(mag(U_rel[2])/(mag(U_plane[1])+SMALL));

			alpha_angle_res[1] = Foam::cos(alpha_angle[1]); //
			beta_angle_res[1] = Foam::cos(beta_angle[1]);
			U_plane[2] = Foam::sqrt(U_rel[0]*U_rel[0] + U_rel[1]*U_rel[1]);
			alpha_angle[2] = SAFE_ACOS(mag(U_rel[0])/(mag(U_plane[2])+SMALL));
			beta_angle[2]  = SAFE_ACOS(mag(U_rel[1])/(mag(U_plane[2])+SMALL));
			alpha_angle_res[2] = Foam::cos(alpha_angle[2]); //
			beta_angle_res[2] = Foam::cos(beta_angle[2]);
			//##############################################################################
			//#########################        WALL FUNCTION     ###########################
			//##############################################################################
			//Pout << "nu[i] = " << nu[i] << endl;
			//Pout << "U_plane = " << U_plane << endl;
			//Pout << "D = " << D << endl;
			for (int q =0; q<3; q++)
			{
				vector D_abs (0,0,0);
				D_abs[0] = mag(D[0]);
				D_abs[1] = mag(D[1]);
				D_abs[2] = mag(D[2]);
				double max_D = std::max(std::max(std::abs(D_abs.x()), std::abs(D_abs.y())), std::abs(D_abs.z()));
				if (D_abs[q] > 0.1*max_D)
				{
					scalar tau_temp = nu[i]*(mag(U_plane[q]) / (mag(D[q])+SMALL));	//**************************
					scalar U_star   = Foam::sqrt(tau_temp/rho[i]);
                    			Yplus[q]   	= mag(D[q])*rho[i]*U_star/nu[i]; 
					scalar ds    	= D50.value();
                   			scalar Kspl    	= rho[i]*ds*U_star/nu[i];
					scalar U_star_new = 0;
					scalar term3 = 0;
					scalar term4 = 0;
					if (Yplus[q] >= 5)
					{
						for (int ii = 0; ii < 40; ++ii) 
                            			{
							scalar term1 = (1 / 0.41) * std::log((9.81 * mag(D[q]) * U_star) / nu[i]);
                                			scalar term2 = (mag((U_plane[q])) / U_star);
                                			if (Kspl>=2.25)
                                			{
                                				term3 = (1 / 0.41) * std::log(1+((0.5 * ds*U_star) / nu[i]));
                                				term4 = ((0.5*ds*rho[i])/(0.5*ds*rho[i]*U_star+rho[i]*nu[i]));
                                			}
                                			else if (Kspl<2.25)
                                			{
                                				term3 = 0;
                                				term4 = 0;
                                			}
                                			scalar numerator = term1 - term2 - term3;
                                			scalar denominator = (1 / (0.41 * U_star)) + (mag(U_plane[q]) / (U_star * U_star)) - term4 + SMALL; //**************************
                                			U_star_new = U_star - (numerator / denominator); 
                                			U_star = (U_star_new);
						
						}
 					U_star_plane[q] = U_star;
                			}
					else
					{
						if (D[q]>0)
						{
						U_star_plane[q] = U_star_plane[q] = rho[i] * Foam::sqrt(nu[i] * mag(U_plane[q]) / D[q]);
						}
					}
				}
				else
				{
					U_star_plane[q] = 0;
				}
			}
			U_star_dir[0] = U_star_dir[0] + U_star_plane[1]*alpha_angle_res[1]+U_star_plane[2]*alpha_angle_res[2];
			U_star_dir[1] = U_star_dir[1] + U_star_plane[0]*alpha_angle_res[0]+U_star_plane[2]*beta_angle_res[2];
			U_star_dir[2] = U_star_dir[2] + U_star_plane[0]*beta_angle_res[0]+U_star_plane[1]*beta_angle_res[1];
		}

	}
	scalar pol = 1;
	for(int hh = 0; hh<3; ++hh)
	{
		if (tau_vis[i][hh] < 0)
		{
			pol = -1;
		}
		tau_vis[i][hh] = pol * U_star_dir[hh] * U_star_dir[hh] * rho[i];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////                          Mean Bed Level                 //////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
scalar sum_Bed_height = 0.0;
vector sum_tau_average = vector::zero;
vector max_tau_max = vector::zero;
vector max_tau_max_each = vector::zero;
scalar nnn = 0;

scalar sum_Bed_height_sq = 0;
scalar dum = 0;
vector sum_tau_sq = vector::zero;
vector dum1 = vector::zero;
vector Theta_sd = vector::zero;
scalar Bed_height_sd = 0;
scalar Bed_height_average = 0;
vector tau_average = vector::zero;
vector Theta_average = vector::zero;
vector Theta_max = vector::zero;
vector Theta_max_each = vector::zero;
const vectorField& cellCentres = mesh.C();
const auto& tau_vis_ref = tau_vis();
forAll(alpha, i)
{
    	if (neeta[i] > 0)
    	{
		nnn = nnn+1;
		sum_Bed_height = sum_Bed_height + cellCentres[i][1];
		sum_tau_average = sum_tau_average + tau_vis_ref[i];
		vector dumm = tau_vis_ref[i];
		for (int ii = 0; ii < 3; ++ii)
		{
			if (dumm[ii]>max_tau_max_each[ii])
			{
				max_tau_max_each[ii] = dumm[ii];
			}
		}
		if (mag(tau_vis_ref[i])>mag(max_tau_max))
		{
			max_tau_max = tau_vis_ref[i];
		}  
	}
}

if (nnn>0)
{
	Bed_height_average = sum_Bed_height/nnn;
	tau_average = (sum_tau_average/nnn);
	forAll(alpha, i)
	{
		if (neeta[i] > 0)
    		{
			dum = Bed_height_average - cellCentres[i][1];
			sum_Bed_height_sq = sum_Bed_height_sq + dum*dum;
			dum1 = tau_average - tau_vis_ref[i];
			sum_tau_sq[0] = sum_tau_sq[0] + dum1[0]*dum1[0];
			sum_tau_sq[1] = sum_tau_sq[1] + dum1[1]*dum1[1];
			sum_tau_sq[2] = sum_tau_sq[2] + dum1[2]*dum1[2];
		}
	}
	Bed_height_sd = Foam::sqrt(max(sum_Bed_height_sq / nnn, scalar(0)));
	Theta_average = tau_average/(mag(g.value())*1.65*D50.value());
	Theta_max = max_tau_max/(mag(g.value())*1.65*D50.value());
	Theta_max_each = max_tau_max_each/(mag(g.value())*1.65*D50.value());
	Theta_sd[0] = Foam::sqrt(max(sum_tau_sq[0] / nnn, scalar(0)))/(mag(g.value())*1.65*D50.value());
	Theta_sd[1] = Foam::sqrt(max(sum_tau_sq[1] / nnn, scalar(0)))/(mag(g.value())*1.65*D50.value());
	Theta_sd[2] = Foam::sqrt(max(sum_tau_sq[2] / nnn, scalar(0)))/(mag(g.value())*1.65*D50.value());
}

Info<< "Bed Level Height mean = "<< Bed_height_average <<endl; 
Info<< "Bed Level Height SD = "<< Bed_height_sd <<endl; 
Info<< "Theta mean = "<< Theta_average <<endl; 
Info<< "Theta max magnitude = "<< Theta_max <<endl; 
Info<< "Theta max (envolope) = "<< Theta_max_each <<endl; 
Info<< "Theta SD = "<< Theta_sd <<endl; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
scalar sum_Bed_height_inlet = 0;
scalar nnnn= 0;
forAll(alpha, i)
{
	if (neeta[i] > 0)
	{
		if (cellCentres[i][0] < 3*D50.value())
		{
			sum_Bed_height_inlet = sum_Bed_height_inlet + cellCentres[i][1];
			nnnn = nnnn + 1;
		}
	}
}
scalar Bed_height_inlet_new = Bed_height_average;
if (nnnn>1)
{
	Bed_height_inlet_new = max(sum_Bed_height_inlet/nnnn,Bed_height_average);
}
Info<< "Bed height inlet before relaxation = "<< Bed_height_inlet_new <<endl; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const scalar relaxFactor = 0.05;

// === Read previous value from registry ===
scalar Bed_height_inlet_old = Bed_height_average;
if (mesh.foundObject<uniformDimensionedScalarField>("Bed_height_inlet"))
{
    Bed_height_inlet_old = mesh.lookupObject
        <uniformDimensionedScalarField>
        ("Bed_height_inlet").value();
}

// === Apply relaxation ===
scalar Bed_height_inlet =
    (1.0 - relaxFactor) * Bed_height_inlet_old
  + relaxFactor         * Bed_height_inlet_new;

Info<< "Bed height inlet after relaxation = "<< Bed_height_inlet<<endl; 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




// Store Bed_height_average in registry for use by BCs
if (mesh.foundObject<uniformDimensionedScalarField>("Bed_height_inlet"))
{
    const_cast<uniformDimensionedScalarField&>(
        mesh.lookupObject<uniformDimensionedScalarField>("Bed_height_inlet")
    ).value() = Bed_height_inlet;
}
else
{
    uniformDimensionedScalarField* bhPtr =
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "Bed_height_inlet",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("Bed_height_inlet", dimLength, Bed_height_inlet)
        );
    bhPtr->store(bhPtr);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////        Store Initial Bed Height ONCE (for mass conservation)   ///////////////////////////////////////////////////////////////////////
///////////////        Uses BOTH registry (fast) and file (restart-robust)     ///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* 
if (!mesh.foundObject<uniformDimensionedScalarField>("Bed_height_inlet_initial"))
{
    // === Store in object registry ===
    // This is the primary mechanism — read by coded BCs each timestep
    uniformDimensionedScalarField* zBed0Ptr =
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "Bed_height_inlet_initial",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar
            (
                "Bed_height_inlet_initial",
                dimLength,
                Bed_height_inlet    // captured at first timestep this runs
            )
        );
    zBed0Ptr->store(zBed0Ptr);
 
    // === Also write to file ===
    // Fallback for restarts — registry is cleared on restart
    // but file persists on disk
    OFstream initFile
    (
        mesh.time().path()/"initialBedHeight.dat"
    );
    initFile << Bed_height_inlet << nl;
 
    Info << "========================================" << nl;
    Info << "Stored Bed_height_inlet_initial = "
         << Bed_height_inlet << " (registry + file)" << nl;
    Info << "========================================" << nl;
}*/