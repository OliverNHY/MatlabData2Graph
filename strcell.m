function cellB=strcell(cellA)
%There is only one data for each element in cellA
%Convert each element to Char class in cellA.
[rs,cs]=size(cellA);
cellB=cell(rs,cs);
breakFlag=1;
indexR=1;
while indexR~=(rs+1) && breakFlag
    indexC=1;
    while indexC~=(cs+1) && breakFlag
        if isa(cellA{indexR,indexC},'double')
            cellB{indexR,indexC}=num2str(cellA{indexR,indexC});
        elseif isa(cellA{indexR,indexC},'cell')
            fprintf('\n\tstrCell: class(cellA{}) should not be cell.');
            breakFlag=0;
        else
            cellB{indexR,indexC}=cellA{indexR,indexC};
        end
        indexC=indexC+1;
    end
    indexR=indexR+1;
end
end