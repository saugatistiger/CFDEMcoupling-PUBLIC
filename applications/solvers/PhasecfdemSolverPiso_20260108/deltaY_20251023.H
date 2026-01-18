neeta = dimensionedScalar("zero", dimless, 0.0);
const volScalarField& nu = mesh.lookupObject<volScalarField>("nu");
forAll ( alpha, i)
{
	if (alpha[i] < 0.73)
	{
    		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			if (alpha[adj[j]]<0.73)
			{
				neeta[i] = 1.0;
				break;
			}
		}
		scalar n = 0;
		for (int j=0; j < adj.size(); ++j)
		{
			if (alpha[adj[j]]<0.73)
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
		//Info << neeta[i] << endl;
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j<adj.size();++j)
		{
			if (alpha[adj[j]] > 0.73)
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

const volScalarField& k = turbulence->k();

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
				double max_D = std::max(std::max(std::abs(D_abs.x()), std::abs(D_abs.y())), std::abs(D_abs.z()));
				if (D_abs[q] > 0.1*max_D)
				{
					scalar tau_temp = nu[i]*(mag(U_plane[q]) / mag(D[q]));	
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
                                			scalar denominator = (1 / (0.41 * U_star)) + (mag(U_plane[q]) / (U_star * U_star)) - term4;
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