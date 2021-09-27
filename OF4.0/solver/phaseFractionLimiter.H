{
/*	
	>> phase fraction limiter is used to limit the disperesed phase fraction between 0 and alphaMax
	>> values exceeding alphaMax are unphysical and cause sharp increase in the particle pressure  
	   and lead to numerical instabilities
*/   
	
    // Calculating exceeding volume fractions
    volScalarField alphaEx(max(alphad - alphaMax+scalar(0.02), scalar(0)));

    // Finding neighbouring cells of the whole domain
    labelListList neighbour = mesh.cellCells(); //2D array
    scalarField cellVolumes(mesh.cellVolumes());

    forAll(alphaEx, celli)
    {
        // Finding the labels of the neighbouring cells for current cell
        labelList neighbourCell = neighbour[celli]; //1D array

        // Initializing neighbouring cells contributions
        scalar neighboursEx = 0.0;	//+ve contribution

        forAll(neighbourCell, cellj)
        {
            labelList neighboursNeighbour = neighbour[neighbourCell[cellj]];
            scalar neighboursNeighbourCellVolumes = 0.0;

            forAll(neighboursNeighbour, cellk)
            {
                neighboursNeighbourCellVolumes +=
                    cellVolumes[neighboursNeighbour[cellk]];
            }

            neighboursEx +=
                alphaEx[neighbourCell[cellj]]*cellVolumes[celli]
               /neighboursNeighbourCellVolumes;
        }

        alphad[celli] += neighboursEx - alphaEx[celli];
    }
	
    alphad.correctBoundaryConditions();

}
