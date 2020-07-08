%Import *.ffe(monoRCS) from dir specified
%Export the 2D origin format dat.
% function: originCell()    cols2mat()
%20200708
clear;clc;tic;
ffeDir='E:\ZM\0Work\3simuModel\SModel\202005825Siks11YZmonoRCS\ffe_116yz\';
resultDir=[ffeDir,'ffe2Origin\'];
if  ~exist(resultDir,'dir')
    mkdir(resultDir)
end
% Exract lines only exist numbers to txt file from ffe file
ffeListOb=dir([ffeDir,'*.ffe']);
ffeFoldnames={ffeListOb.folder};
ffeNames={ffeListOb.name};
ffeLongNames=fullfile(ffeFoldnames,ffeNames);
[ffeFiltFolder,ffeFiltLNames]=filtLine(ffeLongNames,'[#,*,NAN]','filtData');
nTxt=length(ffeFiltLNames);count=0;
for iTxt=1:nTxt
    longName=ffeFiltLNames{iTxt};
    [~,simName,~]=fileparts(longName);
%     f0=simName(strfind(simName,'Fre')+3:strfind(simName,'M')-1);
    obTxt=importdata(longName);
    if isstruct(obTxt)
        data=obTxt.data(:,[2,1,end]);
    else
        data=obTxt(:,[2,1,end]);
    end
    thetaVec=unique(data(:,2));phiVec=unique(data(:,1));
    nr=length(phiVec);nc=length(thetaVec);
%     rcs=NaN(nr,nc);
    [rcs,vecsCell]=cols2mat(data);
    zeroCmpBool=(rcs<0);
    if sum(zeroCmpBool(:))~=0
        rcsdB=rcs;
    else
        rcsdB=20.*log10(rcs);
    end
    varNames={'Phi','RCS'};
    varUnits={'Deg','dBsm'};
    varComnts=num2cell([simName,string([repmat('Theta=',length(vecsCell{2}),1),num2str(vecsCell{2})])']);
    outCell=originCell(varNames,varUnits,varComnts,vecsCell{1},rcsdB);
    writecell(outCell,[resultDir,simName,'.dat']);fprintf('\n%d/%d writecell to\n%s',iTxt,nTxt,[resultDir,simName,'.dat']);
end
elapedTime=toc;
fprintf('\nElaped time: \n(%f s)(%f m or %f h)',elapedTime,elapedTime/60,elapedTime/60/60);
sound(sin(2*pi*25*(1:4000)/500));
