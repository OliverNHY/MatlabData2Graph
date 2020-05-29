clear;clc;
dirName='E:\ZM\0Work\3simuModel\SModel\202005825monoRCS\';
newFolder='format2Org\';
newDir=[dirName,newFolder];
if ~exist(newDir,'dir')
    mkdir(newDir);
end
dirOb=dir([dirName,'*.dat']);
nFiles=length(dirOb);
datNames={dirOb.name};
datDir={dirOb.folder};
longNames=fullfile(datDir,datNames);
for index=1:nFiles
    datLName=longNames{index};
    dataOb=importdata(datLName);
    rawStr=dataOb.colheaders;
    nRaw=length(rawStr);
    varComments=cell(1,nRaw);
    varComments(1)={'RCS'};
    for varIndex=2:nRaw
        if contains(rawStr{varIndex},'Total') || contains(rawStr{varIndex},'Plane')
            varComments(varIndex)=extractBetween(rawStr{varIndex},'- ',' [');
        else
            varComments(varIndex)=extractBetween(rawStr{varIndex},'"',' [');
        end
    end
    varNames(1,1)=extractBetween(rawStr{1},'"Plane Wave ','"');
    varNames(1,2)={'RCS'};
    varUnits={'Deg','dBm^2'};
    xVec=dataOb.data(:,1);data=dataOb.data(:,2:end);
    orgFData=originCell(varNames,varUnits,varComments,xVec,data);
    writecell(orgFData,[newDir,datNames{index}],'Delimiter','tab');
    fprintf('\n%d/%d format %s',index,nFiles,datNames{index});
end
fprintf('\nDone!!!!!!!!!');