% filt string in a text document.
% output the result in a similar name.
function outLongName=filtLine(longName,str)
fpIn=fopen(longName,'r');
[folderName,fileName,ext]=fileparts(longName);
newName=[fileName,'(filt_',str,')',ext];
outLongName=fullfile(folderName,newName);
fpOut=fopen(outLongName,'w');
while ~feof(fpIn)
    lineStr=fgets(fpIn);
    if ~contains(lineStr,str)
        fprintf(fpOut,'%s',lineStr);
    end
end
fclose(fpIn);fclose(fpOut);
end