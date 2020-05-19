% filt string in a text document.
% output the result in a similar name.
function [outDirName,outLongNames]=filtLine(longNames,str,newfolder)
count=0;
outLongNames=cell(length(longNames),1);
for longNameCell=longNames
    count=count+1;
    longName=longNameCell{1};
    fpIn=fopen(longName,'r');
    [folderName,fileName,ext]=fileparts(longName);
    newName=[fileName,ext];%,'(filt_',str,')'
    newDir=fullfile(folderName,newfolder);
    if ~exist(newDir,'dir')
        mkdir(newDir);
    end
    outLongName=fullfile(newDir,newName);
    outLongNames{count}=outLongName;
    fpOut=fopen(outLongName,'w');
    while ~feof(fpIn)
        lineStr=fgets(fpIn);
        if ~contains(lineStr,str)
            fprintf(fpOut,'%s',lineStr);
        end
    end
    fclose(fpIn);fclose(fpOut);
end
outDirName=newDir;
fprintf('\nfiltLine(): Filt %s for %d files.\n',str,count);
end