function writeCellTo(dir,name,cellA)
%%There is only one data for each element in cellA
[rs,cs]=size(cellA);
breakFlag=1;
fileID=fopen([dir,name],'w');
if fileID==-1
    breakFlag=0;
    fprintf('\n\twriteCellTo(): Failed to open %s\n',[dir,name]);
end
indexR=1;
while indexR~=1+rs && breakFlag
    indexC=1;
    while indexC~=1+cs && breakFlag
        if isa(cellA{indexR,indexC},'double')
            fprintf(fileID,'%f ',cellA{indexR,indexC});
        elseif isa(cellA{indexR,indexC},'char')
            fprintf(fileID,'%s ',cellA{indexR,indexC});
        else
            fprintf('\n\twriteCellTo(): class(cellA{}) is neither Double nor Char (May be Cell)\n');
            breakFlag=0;
        end
        indexC=indexC+1;
    end
    fprintf(fileID,'\n');
    indexR=indexR+1;
end
if fileID~=-1
    fclose(fileID);
    fprintf('\n\twriteCell(): Write to %s\n',[dir,name]);
end
end
            
        
