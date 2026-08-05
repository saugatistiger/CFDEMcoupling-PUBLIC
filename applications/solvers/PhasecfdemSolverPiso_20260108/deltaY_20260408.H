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
			if (alpha[adj[j]]<0.76)
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
forAll ( alpha, i)
{
	if (((neeta[i] > 0) && (mag(Us[i])>4*u_star.value())))
	{
		neeta[i] = 0;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
forAll ( alpha, i)
{
	if (neeta[i] > 0)
	{
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			vector DDD = mesh.C()[adj[j]] - mesh.C()[i];
			if ((mag(DDD[0])>0.001) && (DDD[1]<0.001) && (mag(DDD[2])<0.001))
			{
				neeta[adj[j]] = 1;
			}
		}
	}
}

forAll ( alpha, i)
{
	if (neeta[i] > 0)
	{
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			vector DDD = mesh.C()[adj[j]] - mesh.C()[i];
			if ((mag(DDD[0])>0.001) && (DDD[1]<0.001) && (mag(DDD[2])<0.001))
			{
				neeta[adj[j]] = 1;
			}
		}
	}
}
forAll ( alpha, i)
{
	if (neeta[i] > 0)
	{
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			vector DDD = mesh.C()[adj[j]] - mesh.C()[i];
			if ((mag(DDD[0])>0.001) && (DDD[1]<0.001) && (mag(DDD[2])<0.001))
			{
				neeta[adj[j]] = 1;
			}
		}
	}
}
forAll ( alpha, i)
{
	if (neeta[i] > 0)
	{
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			vector DDD = mesh.C()[adj[j]] - mesh.C()[i];
			if ((mag(DDD[0])>0.001) && (DDD[1]<0.001) && (mag(DDD[2])<0.001))
			{
				neeta[adj[j]] = 1;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
forAll ( alpha, i)
{
	if (neeta[i] > 0)
	{
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			vector DDD = mesh.C()[adj[j]] - mesh.C()[i];
			if ((mag(DDD[0])<0.001) && (DDD[1]>0.001) && (mag(DDD[2])<0.001))
			{
				neeta[adj[j]] = 0;
				//////////////////////////////////////////////////////////////////////////////
				labelList adj_j = mesh.cellCells()[adj[j]];
				for (int jj = 0; jj < adj_j.size(); ++jj)
				{
					vector DDDD = mesh.C()[adj_j[jj]] - mesh.C()[adj[j]];
					if ((mag(DDDD[0])<0.001) && (DDDD[1]>0.001) && (mag(DDDD[2])<0.001))
					{
						neeta[adj_j[jj]] = 0;
						//////////////////////////////////////////////////////////////
						labelList adj_jj = mesh.cellCells()[adj_j[jj]];
						for (int jjj = 0; jjj < adj_jj.size(); ++jjj)
						{
							vector DDDDD = mesh.C()[adj_jj[jjj]] - mesh.C()[adj_j[jj]];
							if ((mag(DDDDD[0])<0.001) && (DDDDD[1]>0.001) && (mag(DDDDD[2])<0.001))
							{
								neeta[adj_jj[jjj]] = 0;
								////////////////////////////////////////////////////////////////////////////////////
								labelList adj_jjj = mesh.cellCells()[adj_jj[jjj]];
								for (int jjjj = 0; jjjj < adj_jjj.size(); ++jjjj)
								{
									vector DDDDDD = mesh.C()[adj_jjj[jjjj]] - mesh.C()[adj_jj[jjj]];
									if ((mag(DDDDDD[0])<0.001) && (DDDDDD[1]>0.001) && (mag(DDDDDD[2])<0.001))
									{
										neeta[adj_jjj[jjjj]] = 0;
										////////////////////////////////////////////////////////////////////////////////////
										labelList adj_jjjj = mesh.cellCells()[adj_jjj[jjjj]];
										{
											for (int jjjjj = 0; jjjjj < adj_jjjj.size(); ++jjjjj)
											{
												vector DDDDDDD = mesh.C()[adj_jjjj[jjjjj]] - mesh.C()[adj_jjj[jjjj]];
												if ((mag(DDDDDDD[0])<0.001) && (DDDDDDD[1]>0.001) && (mag(DDDDDDD[2])<0.001))
												{
													neeta[adj_jjjj[jjjjj]] = 0;
													////////////////////////////////////////////////////////////////////////////////////
													labelList adj_jjjjj = mesh.cellCells()[adj_jjjj[jjjjj]];
													{
														for (int jjjjjj = 0; jjjjjj < adj_jjjjj.size(); ++jjjjjj)
														{
															vector DDDDDDDD = mesh.C()[adj_jjjjj[jjjjjj]] - mesh.C()[adj_jjjj[jjjjj]];
															if ((mag(DDDDDDDD[0])<0.001) && (DDDDDDDD[1]>0.001) && (mag(DDDDDDDD[2])<0.001))
															{
																neeta[adj_jjjjj[jjjjjj]] = 0;
																////////////////////////////////////////////////////////////////////////////////////
																labelList adj_jjjjjj = mesh.cellCells()[adj_jjjjj[jjjjjj]];
																{
																	for (int jjjjjjj = 0; jjjjjjj < adj_jjjjjj.size(); ++jjjjjjj)
																	{
																		vector DDDDDDDDD = mesh.C()[adj_jjjjjj[jjjjjjj]] - mesh.C()[adj_jjjjj[jjjjjj]];
																		if ((mag(DDDDDDDDD[0])<0.001) && (DDDDDDDDD[1]>0.001) && (mag(DDDDDDDDD[2])<0.001))
																		{
																			neeta[adj_jjjjjj[jjjjjjj]] = 0;
																			////////////////////////////////////////////////////////////////////////////////////
																			labelList adj_jjjjjjj = mesh.cellCells()[adj_jjjjjj[jjjjjjj]];
																			{
																				for (int jjjjjjjj = 0; jjjjjjjj < adj_jjjjjjj.size(); ++jjjjjjjj)
																				{
																					vector DDDDDDDDDD = mesh.C()[adj_jjjjjjj[jjjjjjjj]] - mesh.C()[adj_jjjjjj[jjjjjjj]];
																					if ((mag(DDDDDDDDDD[0])<0.001) && (DDDDDDDDDD[1]>0.001) && (mag(DDDDDDDDDD[2])<0.001))
																					{
																						neeta[adj_jjjjjjj[jjjjjjjj]] = 0;
																					}
																				}
																			}
																			////////////////////////////////////////////////////////////////////////////////////
																		}
																	}
																}
																////////////////////////////////////////////////////////////////////////////////////
															}
														}
													}
													////////////////////////////////////////////////////////////////////////////////////
												}
											}
										}
										////////////////////////////////////////////////////////////////////////////////////
									}
								}
								/////////////////////////////////////////////////////////////////////////////////////
							}
						}
						//////////////////////////////////////////////////////////////
					}	
				}
				//////////////////////////////////////////////////////////////////////////////
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

forAll ( alpha, i)
{
	if (neeta[i] > 0)
	{
		scalar nnn = 0;
		labelList adj = mesh.cellCells()[i];
		for (int j = 0; j < adj.size(); ++j)
		{
			if (neeta[adj[j]] > 0)
			{
				nnn = nnn + 1;
			}
		}
		if (nnn==0)
		{
			neeta[i] = 0;
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
		//Info << neeta[i] << endl;
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
			U_plane[0] = Foam::sqrt(U_rel[1]*U_rel[1] + U_rel[2]*U_rel[2]); //yz-plane
			alpha_angle[0] = Foam::acos(mag(U_rel[1])/(mag(U_plane[0])+SMALL)); //**************************
			beta_angle[0] = Foam::acos(mag(U_rel[2])/(mag(U_plane[0])+SMALL)); // **************************
			alpha_angle_res[0] = Foam::cos(alpha_angle[0]); //
			beta_angle_res[0] = Foam::cos(beta_angle[0]);
			U_plane[1] = Foam::sqrt(U_rel[0]*U_rel[0] + U_rel[2]*U_rel[2]); //xz-plane
			alpha_angle[1] = Foam::acos(mag(U_rel[0])/(mag(U_plane[1])+SMALL)); //**************************
			beta_angle[1] = Foam::acos(mag(U_rel[2])/(mag(U_plane[1])+SMALL));  //**************************
			alpha_angle_res[1] = Foam::cos(alpha_angle[1]); //
			beta_angle_res[1] = Foam::cos(beta_angle[1]);
			U_plane[2] = Foam::sqrt(U_rel[0]*U_rel[0] + U_rel[1]*U_rel[1]); //yz-plane
			alpha_angle[2] = Foam::acos(mag(U_rel[0])/mag(U_plane[2])+SMALL);  //**************************
			beta_angle[2] = Foam::acos(mag(U_rel[1])/mag(U_plane[2])+SMALL);   //**************************
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


// 1. Before the block - check fields are valid
Info<< "=== Bed Stats START ===" << endl;
Info<< "tau_vis[0] = " << tau_vis_ref[0] << endl;
Info<< "cellCentres[0] = " << cellCentres[0] << endl;
Info<< "neeta[0] = " << neeta[0] << endl;


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


// 2. After first forAll - check accumulation is sane
Info<< "nnn = " << nnn << endl;
Info<< "sum_Bed_height = " << sum_Bed_height << endl;
Info<< "sum_tau_average = " << sum_tau_average << endl;
Info<< "max_tau_max = " << max_tau_max << endl;


if (nnn>0)
{
	Bed_height_average = sum_Bed_height/nnn;
	tau_average = (sum_tau_average/nnn);

	
	// 3. After tau_average is computed - check for NaN/Inf
	Info<< "Bed_height_average = " << Bed_height_average << endl;
	Info<< "tau_average = " << tau_average << endl;


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


// 4. After second forAll - check sums before sqrt
Info<< "sum_Bed_height_sq = " << sum_Bed_height_sq << endl;
Info<< "sum_tau_sq = " << sum_tau_sq << endl;



Info<< "Bed Level Height mean = "<< Bed_height_average <<endl; 
Info<< "Bed Level Height SD = "<< Bed_height_sd <<endl; 
Info<< "Theta mean = "<< Theta_average <<endl; 
Info<< "Theta max magnitude = "<< Theta_max <<endl; 
Info<< "Theta max (envolope) = "<< Theta_max_each <<endl; 
Info<< "Theta SD = "<< Theta_sd <<endl; 


// 5. After the block - confirm solver fields still intact
Info<< "=== Bed Stats END ===" << endl;
Info<< "tau_vis_ref[0] post = " << tau_vis_ref[0] << endl;
Info<< "cellCentres[0] post = " << cellCentres[0] << endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

