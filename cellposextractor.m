function [cellpositions celltype]= cellposextractor(psim)

    number_of_cells = psim.ReturnNumberCells;
       
    for index = 1:number_of_cells
        cellpositions(index,1) = psim.ReturnCellPositions(index*2-1-1);
        cellpositions(index,2) = psim.ReturnCellPositions(index*2-1);
        celltype(index) = psim.ReturnCellType(index-1);
    end

end