// code found from https://www.cfd-online.com/Forums/openfoam-programming-development/80038-cell-size-x-y-z.html
//	const surfaceVectorField& deltas=mesh.delta();
//	const volVectorField& del=mesh.delta();
//	Info<<deltas.component(0)<<endl;
//	Info<<deltas<<endl;
//////////////////////////////////////////////////////////////////////////////////////////
//std::vector<scalar> xDims(mesh.C().size(), 0);
//volScalarField xDims("p",0.*p);
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
