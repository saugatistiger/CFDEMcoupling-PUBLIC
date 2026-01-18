neeta = dimensionedScalar("zero", dimless, 0.0);
const volScalarField& nu = mesh.lookupObject<volScalarField>("nu");
forAll ( alpha, i)
{
	if (alpha[i] < 1)
	{
    		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			if (alpha[adj[j]]>0.73)
			{
				neeta[i] = 1.0;
				break;
			}
		}
		scalar n = 0;
		for (int j=0; j < adj.size(); ++j)
		{
			if (alpha[adj[j]]>0.73)
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

labelList 	innerList(6,-1);
labelListList 	neeta_nei(alpha.size(),innerList); // neighbour of the neeta

forAll ( alpha, i)
{
	int n = 0;
	if ( neeta[i] == 1 )
	{
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j<adj.size();++j)
		{
			if (alpha[adj[j]] > 0.73)
			{
				neeta_nei[i][n]=adj[j];
				n = n +1;
			}
		}
	}
}
Pout << "alpha =" <<alpha.size() << endl;

forAll (alpha, i)
{
	vector U_star_dir(0,0,0);
	if ( neeta[i] == 1 )
	{
	int a = 0;
	for (int k = 0; k < 6 ; ++k)
	{
		if (neeta_nei[i][k]==-1)
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
		for (int j=0 ; j < a ; ++j)
		{
			D = mesh.C()[neeta_nei[i][j]] - mesh.C()[i];
			U_rel = U[neeta_nei[i][j]] - U[i] ;
			U_plane[0] = Foam::sqrt(U_rel[1]*U_rel[1] + U_rel[2]*U_rel[2]); //yz-plane
			alpha_angle[0] = Foam::acos(U_rel[1]/U_plane[0]); //
			beta_angle[0] = Foam::acos(U_rel[2]/U_plane[0]);
			alpha_angle_res[0] = Foam::cos(alpha_angle[0]); //
			beta_angle_res[0] = Foam::cos(beta_angle[0]);
			U_plane[1] = Foam::sqrt(U_rel[0]*U_rel[0] + U_rel[2]*U_rel[2]); //xz-plane
			alpha_angle[1] = Foam::acos(U_rel[0]/U_plane[1]);
			beta_angle[1] = Foam::acos(U_rel[2]/U_plane[1]);
			alpha_angle_res[1] = Foam::cos(alpha_angle[1]); //
			beta_angle_res[1] = Foam::cos(beta_angle[1]);
			U_plane[2] = Foam::sqrt(U_rel[0]*U_rel[0] + U_rel[1]*U_rel[1]); //yz-plane
			alpha_angle[2] = Foam::acos(U_rel[0]/U_plane[2]);
			beta_angle[2] = Foam::acos(U_rel[1]/U_plane[2]);
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
				double max_D = *std::max_element(D_abs.begin(), D_abs.end());
				if (D_abs[q] > 0.1*max_D)
				{
					scalar tau_temp = nu[i]*(mag(U_plane[q]) / mag(D[q]));	
					scalar U_star   = Foam::sqrt(tau_temp/rho[i]);
                        		Yplus[q]   	= mag(D[q])*rho[i]*U_star/nu[i]; 
					scalar ds    	= 0.01;
                        		scalar Kspl    	= rho[i]*2.5*ds*U_star/nu[i];
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
                                				term3 = (1 / 0.41) * std::log(1+((0.5 * ds*2.5 * U_star) / nu[i]));
                                				term4 = ((0.5*ds*2.5*rho[i])/(0.5*2.5*ds*rho[i]*U_star+rho[i]*nu[i]));
                                			}
                                			else if (Kspl<2.25)
                                			{
                                				term3 = 0;
                                				term4 = 0;
                                			}
                                			scalar numerator = term1 - term2 - term3;
                                			scalar denominator = (1 / (0.41 * U_star)) + (mag(U_plane[q]) / (U_star * U_star)) - term4;
                                			U_star_new = U_star - (numerator / denominator); 
                                			U_star = (U_star_new);
						
						}
 					U_star_plane[q] = U_star;
                			}
					else
					{
						U_star_plane[q] = U_star_plane[q] = rho[i] * Foam::sqrt(nu[i] * mag(U_plane[q]) / D[q]);
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
//	Info<<xDims<<endl;
}
//Info<<xDims<<endl;
