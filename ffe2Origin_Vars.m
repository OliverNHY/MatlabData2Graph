% Get N-D Array
% Import *.ffe(monoRCS) from dir specified
% Extract key data and save as Origin data format.
% Function setLinWidth()   plot3D()    filtLine()
% 20200526
clear;clc;
tic;
varStrLimit={'L','H';'H','A1';'A1','A2';'A2',')(Fre';'Fre','M';'Pol',')'};[nVar,~]=size(varStrLimit);
ffeDir='E:\ZM\0Work\3simuModel\20200523GaussModel\';
matName='data.mat';matLongName=[ffeDir,matName];
% Filter lines contain characters specified out in ffe;
% result in .\filtData\*.txt
ffeListOb=dir([ffeDir,'*.ffe']);
ffeFoldnames={ffeListOb.folder};
ffeNames={ffeListOb.name};
ffeLongNames=fullfile(ffeFoldnames,ffeNames);
[ffeFiltFolder,ffeFiltLNames]=filtLine(ffeLongNames,'#','filtData');
% Get variate names and its value from (ffe) file names
varsCache=NaN(1,nVar);varsMat=NaN(length(ffeNames),nVar);
phiVec=[];thetaVec=[];
for indexffe=1:length(ffeNames)
    for indexVar=1:nVar
        varsCache(indexVar)=str2double(extractBetween(ffeNames{indexffe},varStrLimit{indexVar,1},varStrLimit{indexVar,2}));
    end
    varsMat(indexffe,:)=varsCache;
    phiStr=extractBetween(ffeNames{indexffe},'_phi','_theta');
    phiVecS=str2double(split(phiStr,{'to','dphi'}));phiVecCache=phiVecS(1):phiVecS(3):phiVecS(2);
    phiVec=unique([phiVec,phiVecCache],'sorted');
    thetaStr=extractBetween(ffeNames{indexffe},'_theta','_Pol');
    thetaVecS=str2double(split(thetaStr,{'to','dtheta'}));thetaVecCache=thetaVecS(1):thetaVecS(3):thetaVecS(2);
    thetaVec=unique([thetaVec,thetaVecCache],'sorted');
end
varsVec=cell(1,nVar);
for indexVar=1:nVar
    varsVec{indexVar}=unique(varsMat(:,indexVar),'sorted');
end
varsVecOb=cell2struct(varsVec,varStrLimit(:,1),2);
varsListOb=cell2struct(num2cell(varsMat),varStrLimit(:,1),2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Extract all variates to one n-dimision array
nTxt=length(ffeFiltLNames);
for count=1:nTxt
    longName=ffeFiltLNames{count};
    [~,simName,~]=fileparts(longName);
    f0=extractBetween(simName,'Fre','M');f0=f0{1};
    obTxt=importdata(longName);
    data=obTxt.data(:,[1,2,end]);
    thetaVecThis=unique(data(:,1));phiVecThis=unique(data(:,2));
    [~,thetaCs]=find(thetaVec==thetaVecThis);%thetaVec=(1,n),thetaVecThis=(m,1)
    [~,phiCs]=find(phiVec==phiVecThis);
    nT=length(thetaVecThis);nP=length(phiVecThis);
    rcs=reshape(data(:,3),[nT,nP]);% Different theta for each same column, different phi for each same row.
    rcsDB=10.*log(rcs)./log(10);
    extDimIndexs=NaN(1,nVar);
    for iVar=1:nVar
        extDimIndexs(1,iVar)=find(varsVec{iVar}==varsMat(count,iVar));
    end
    rcs_ND(thetaCs,phiCs,extDimIndexs(1),extDimIndexs(2),extDimIndexs(3),extDimIndexs(4),extDimIndexs(5),extDimIndexs(6))=rcs;
    rcsDB_ND(thetaCs,phiCs,extDimIndexs(1),extDimIndexs(2),extDimIndexs(3),extDimIndexs(4),extDimIndexs(5),extDimIndexs(6))=rcsDB;
    fprintf('\n%d/%d Absorb %s',count,nTxt,simName);
end
save(matLongName);
fprintf('\nDone!!!!!!!!!!!\n')
toc

%% Select and Plot RCS(L,F).
rcsUnit='dBsm';
% rcsUnit='m^2';
matLongName='E:\ZM\0Work\3simuModel\20200523GaussModel\data.mat';
tic
if exist(matLongName,'file')
    load(matLongName);
end
fontName='Times New Roman';fontWeight='bold';
ExportDir=[ffeDir,'ffe2Origin\'];
if ~exist(ExportDir,'dir')
    mkdir(ExportDir)
end
theta90Index=find(thetaVec==90);
phi0Index=find(phiVec==0);
for iPol=1:length(varsVecOb.Pol) 
    count=0;picsHd=figure;
    set(gcf,'outerposition',get(0,'screensize'));
    titleNamePics=['RCS(L,Fre)','(',rcsUnit,')',num2str(length(varsVecOb.A2)),'Pol',num2str(varsVecOb.Pol(iPol))];
    for iA2=1:length(varsVecOb.A2)
%         clf;
        count=count+1;
        if strcmp(rcsUnit,'dBsm')
            dataN=rcsDB_ND(theta90Index,phi0Index,:,1,1,iA2,:,iPol);
        else
            dataN=rcs_ND(theta90Index,phi0Index,:,1,1,iA2,:,iPol);
        end
        data=reshape(dataN,length(varsVecOb.L),length(varsVecOb.Fre));
        [LGrid,FreGrid]=ndgrid(varsVecOb.L,varsVecOb.Fre);%       +1000 
        np=ceil(sqrt(length(varsVecOb.A2)));
        subplot(np,np,count);
        [picHdDB,axHdDB,cbHdDB]=plot3D(LGrid,FreGrid,data,data);drawnow;
        titleName=['A2=',num2str(abs(varsVecOb.A2(iA2))),'Deg,Pol=',num2str(varsVecOb.Pol(iPol))];
        title(titleName,'fontName',fontName,'fontWeight',fontWeight);
        xlabel('L(mm)','fontName',fontName,'fontWeight',fontWeight);
        ylabel('Frq(MHz)','fontName',fontName,'fontWeight',fontWeight);
        xticks(linspace(min(varsVecOb.L),max(varsVecOb.L),5));yticks(varsVecOb.Fre);
%         xticklabels(num2cell(varsVecOb.L));yticklabels(num2cell(varsVecOb.Fre));
        set(axHdDB,'fontName',fontName,'fontWeight',fontWeight);
        cbHdDB.Label.String=rcsUnit;set(cbHdDB,'fontName',fontName,'fontWeight',fontWeight);
        if varsVecOb.Pol(iPol)==0
            if strcmp(rcsUnit,'dBsm')
                caxis([-45,10]);
            else
                caxis([0,12]);
            end
        else
            if strcmp(rcsUnit,'dBsm')
                caxis([-15,10]);
            else
                caxis([0,6]);
            end
        end
        % hold on;contour(LGrid,FreGrid,data,'black');
        freInterp=linspace(min(varsVecOb.Fre),max(varsVecOb.Fre),100);
        lamdaVec=300./freInterp.*1e3; 
%         hold on;plot(lamdaVec.*4,freInterp,'black','LineWidth',1,'LineStyle','--');
%         hold on;plot(lamdaVec.*2,freInterp,'black','LineWidth',1,'LineStyle','--');
        hold on;plot(lamdaVec,freInterp,'black','LineWidth',1,'LineStyle','--');drawnow;
        hold on;plot(lamdaVec./2,freInterp,'black','LineWidth',1,'LineStyle','--');drawnow;
%         hold on;pltHd=plot(lamdaVec./4,freInterp,'black','LineWidth',1,'LineStyle','--');
        % pltHd.LineWidth=2;
        % set(pltHd,'LineWidth',2)
%         picLongName=[ExportDir,titleName,'.bmp'];
%         saveas(picHdDB,picLongName);
    end
    picsLongName=[ExportDir,titleNamePics,'.bmp'];
    saveas(picsHd,picsLongName);
end
toc

%% Select and Plot RCS(L,A2).
rcsUnit='dBsm';
% rcsUnit='m^2';
matLongName='E:\ZM\0Work\3simuModel\20200523GaussModel\data.mat';
tic
% if exist(matLongName,'file')
%     load(matLongName);
% end
fontName='Times New Roman';fontWeight='bold';
ExportDir=[ffeDir,'ffe2Origin\'];
if ~exist(ExportDir,'dir')
    mkdir(ExportDir)
end
theta90Index=find(thetaVec==90);
phi0Index=find(phiVec==0);
for iPol=1:length(varsVecOb.Pol) 
    count=0;picsHd=figure;
    set(gcf,'outerposition',get(0,'screensize'));
    titleNamePics=['RCS(L,A2)','(',rcsUnit,')',num2str(length(varsVecOb.Fre)),'Pol',num2str(varsVecOb.Pol(iPol))];    
    for iFre=1:length(varsVecOb.Fre)
%         clf;
        count=count+1;
        if strcmp(rcsUnit,'dBsm')
            dataN=rcsDB_ND(theta90Index,phi0Index,:,1,1,:,iFre,iPol);
        else
            dataN=rcs_ND(theta90Index,phi0Index,:,1,1,:,iFre,iPol);
        end
        data=reshape(dataN,length(varsVecOb.L),length(varsVecOb.A2));
        [LGrid,A2Grid]=ndgrid(varsVecOb.L,abs(varsVecOb.A2));%+1000
        np=ceil(sqrt(length(varsVecOb.Fre)));
        subplot(np,np,count);
        [picHdDB,axHdDB,cbHdDB]=plot3D(LGrid,A2Grid,data,data);drawnow;
        titleName=['Fre=',num2str(abs(varsVecOb.Fre(iFre))),'MHz,Pol=',num2str(varsVecOb.Pol(iPol))];
        title(titleName,'fontName',fontName,'fontWeight',fontWeight);
        xlabel('L(mm)','fontName',fontName,'fontWeight',fontWeight);
        ylabel('A2(Deg)','fontName',fontName,'fontWeight',fontWeight);
        xticks(linspace(min(varsVecOb.L),max(varsVecOb.L),5));A2Vec=abs(varsVecOb.A2);yticks(A2Vec(end:-1:1));
        set(axHdDB,'fontName',fontName,'fontWeight',fontWeight);
        cbHdDB.Label.String=rcsUnit;set(cbHdDB,'fontName',fontName,'fontWeight',fontWeight);
        if varsVecOb.Pol(iPol)==0
            if strcmp(rcsUnit,'dBsm')
                caxis([-45,10]);
            else
                caxis([0,12]);
            end
        else
            if strcmp(rcsUnit,'dBsm')
                caxis([-15,10]);
            else
                caxis([0,6]);
            end
        end
        % hold on;contour(LGrid,FreGrid,data,'black');
%         picLongName=[ExportDir,titleName,'.bmp'];
%         saveas(picHdDB,picLongName);
        % frame=getframe(picHdDB);
        % imwrite(frame.cdata,[ExportDir,['Fre=',num2str(varsVecOb.Fre(iFre)),'Pol=',num2str(varsVecOb.Pol(iPol))],'.bmp']);
    end
    picsLongName=[ExportDir,titleNamePics,'.bmp'];
    saveas(picsHd,picsLongName);
end
toc

%% Select and Plot Difference A2.
tic
% matLongName='E:\ZM\0Work\3simuModel\20200523GaussModel\data.mat';
% if exist(matLongName,'file')
%     load(matLongName);
% end
fontName='Times New Roman';fontWeight='bold';
ExportDir=[ffeDir,'ffe2Origin\'];
if ~exist(ExportDir,'dir')
    mkdir(ExportDir)
end
theta90Index=find(thetaVec==90);
phi0Index=find(phiVec==0);
for iPol=1:length(varsVecOb.Pol) 
    count=0;
    for iA2=1:length(varsVecOb.A2)
        clf;
        count=count+1;
        dataN=rcsDB_ND(theta90Index,phi0Index,:,1,1,iA2,:,iPol);
        data=reshape(dataN,length(varsVecOb.L),length(varsVecOb.Fre));
        if count==1
            Ddata=data;
            dataCache=data;
        else
            Ddata=data-dataCache;
        end
        [LGrid,FreGrid]=ndgrid(varsVecOb.L,varsVecOb.Fre);
        [picHdDB,axHdDB,cbHdDB]=plot3D(LGrid,FreGrid,Ddata,double(Ddata>0));
        titleName=['Difference(A45)_-A2=',num2str(abs(varsVecOb.A2(iA2))),'Deg,Pol=',num2str(varsVecOb.Pol(iPol))];
        title(titleName,'fontName',fontName,'fontWeight',fontWeight);
        xlabel('L(mm)','fontName',fontName,'fontWeight',fontWeight);
        ylabel('Frq(MHz)','fontName',fontName,'fontWeight',fontWeight);
        set(axHdDB,'fontName',fontName,'fontWeight',fontWeight);
        cbHdDB.Label.String='dBsm';set(cbHdDB,'fontName',fontName,'fontWeight',fontWeight);
%         if varsVecOb.Pol(iPol)==0
%             caxis([-40,10]);
%         else
%             caxis([-15,10]);
%         end
        % hold on;contour(LGrid,FreGrid,data,'black');
        freInterp=linspace(min(varsVecOb.Fre),max(varsVecOb.Fre),100);
        lamdaVec=300./freInterp.*1e3; 
        hold on;plot(lamdaVec.*4,freInterp,'black','LineWidth',1,'LineStyle','--');
        hold on;plot(lamdaVec.*2,freInterp,'black','LineWidth',1,'LineStyle','--');
        hold on;plot(lamdaVec,freInterp,'black','LineWidth',1,'LineStyle','--');
        hold on;plot(lamdaVec./2,freInterp,'black','LineWidth',1,'LineStyle','--');
        hold on;pltHd=plot(lamdaVec./4,freInterp,'black','LineWidth',1,'LineStyle','--');
        % pltHd.LineWidth=2;
        % set(pltHd,'LineWidth',2)
        picLongName=[ExportDir,titleName,'.bmp'];
        saveas(picHdDB,picLongName);
        % frame=getframe(picHdDB);
        % imwrite(frame.cdata,[ExportDir,['A2=',num2str(varsVecOb.A2(iA2)),'Pol=',num2str(varsVecOb.Pol(iPol))],'.bmp']);
    end
end
toc

%% Select and Plot RCS(L)2D.
rcsUnit='dBsm';
% rcsUnit='m^2';
matLongName='E:\ZM\0Work\3simuModel\20200523GaussModel\data.mat';
tic
% if exist(matLongName,'file')
%     load(matLongName);
% end
fontName='Times New Roman';fontWeight='bold';
ExportDir=[ffeDir,'ffe2Origin\'];
if ~exist(ExportDir,'dir')
    mkdir(ExportDir)
end
theta90Index=find(thetaVec==90);
phi0Index=find(phiVec==0);
A2_v=1;
A20Index=find(abs(varsVecOb.A2)==A2_v);
for iPol=1:length(varsVecOb.Pol) 
    count=0;picsHd=figure;
    titleNamePics=['2D_RCS(L,Fre)','(',rcsUnit,')','A2',num2str(A2_v),'Pol',num2str(varsVecOb.Pol(iPol))];
    for iFre=1:length(varsVecOb.Fre)
%         clf;
        count=count+1;
        if strcmp(rcsUnit,'dBsm')
            dataN=rcsDB_ND(theta90Index,phi0Index,:,1,1,A20Index,iFre,iPol);
        else
            dataN=rcs_ND(theta90Index,phi0Index,:,1,1,A20Index,iFre,iPol);
        end
        data=reshape(dataN,length(varsVecOb.L),1);
        hold on;pltHd=plot(varsVecOb.L,data);
        pltHd.LineWidth=2;
    end
    legendStr=cellfun(@strcat,cellstr(repmat('Fre=',size(varsVecOb.Fre))),cellstr(num2str(varsVecOb.Fre,'%-.1f')),cellstr(repmat('MHz',size(varsVecOb.Fre))), 'UniformOutput',false)';
    legend(legendStr);
    set(pltHd,'LineWidth',2);grid on;
    titleName=['A2=',num2str(A2_v),'Deg,Pol=',num2str(varsVecOb.Pol(iPol))];
    title(titleName,'fontName',fontName,'fontWeight',fontWeight);
    xlabel('L(mm)','fontName',fontName,'fontWeight',fontWeight);
    ylabel(['RCS','(',rcsUnit,')'],'fontName',fontName,'fontWeight',fontWeight);
    xticks(linspace(min(varsVecOb.L),max(varsVecOb.L),9));
    xticklabels(num2cell(linspace(min(varsVecOb.L),max(varsVecOb.L),9)));
    xlim([min(varsVecOb.L),max(varsVecOb.L)]);
    set(gca,'fontName',fontName,'fontWeight',fontWeight);
    picsLongName=[ExportDir,titleNamePics,'.bmp'];
    set(gcf,'outerposition',get(0,'screensize'));
    saveas(picsHd,picsLongName);
end
toc

%% 
%     gcf.OuterPosition=get(0,'screensize');
set(gcf,'outerposition',get(0,'screensize'));
%     figHd=gcf;
frame=getframe(figHd);
imwrite(frame.cdata,[ExportDir,simName,'(2D_3D).bmp']);
vThetaMat=zeros(nT,2*nP);
vThetaMat(:,1:2:end-1)=rcs;vThetaMat(:,2:2:end)=rcsDB;%independent variate is theta
vPhiMat=zeros(nP,2*nT);
vPhiMat(:,1:2:end-1)=rcs';vPhiMat(:,2:2:end)=rcsDB';%independent variate is phi
%%% Data header
%   LongName	Theta           RCS     RCS         ...
%   Units	    Deg             m^2     dBsm        ...
%   Comments	Frquency=1GHz	Phi=0	Phi=0       ...
elem2ThetaVec=repelem(thetaVec,2);elem2PhiVec=repelem(phiVec,2);
theta2Str=cellstr(num2str(elem2ThetaVec,'%-.1f'));phi2Str=cellstr(num2str(elem2PhiVec,'%-.1f'));
vThetaComments=cellfun(@strcat,cellstr(repmat('Phi=',size(phi2Str))),phi2Str, 'UniformOutput',false)';
vPhiComments=cellfun(@strcat,cellstr(repmat('Theta=',size(theta2Str))),theta2Str, 'UniformOutput',false)';
%%%Produce 2D 3D data for Origin
vTheta2D=originCell({'Theta','RCS'},{'Deg','m^2','dBsm'},[['Frequency=',f0,'MHz'],vThetaComments],thetaVec,vThetaMat);
vPhi2D=originCell({'Phi','RCS'},{'Deg','m^2','dBsm'},[['Frequency=',f0,'MHz'],vPhiComments],phiVec,vPhiMat);
rcs3D=num2cell([0,phiVec';thetaVec,rcs]);
rcs3D(1,1)={'Theta\Phi'};
rcsDB3D=num2cell([0,phiVec';thetaVec,rcsDB]);
rcsDB3D(1,1)={'Theta\Phi'};
%%%Export 2D 3D data for Origin
writeCellTo(ExportDir,[simName,'(vTheta2D).txt'],vTheta2D);
writeCellTo(ExportDir,[simName,'(vPhi2D).txt'],vPhi2D);
writeCellTo(ExportDir,[simName,'(dB3D).txt'],rcsDB3D);
writeCellTo(ExportDir,[simName,'(3D).txt'],rcs3D);

%%%Plot 2D and 3D graph
clf;
figHd=gcf;
figHd.Name=simName;
subplot(2,2,1);
vThetaHd=plot(thetaVec,rcsDB);%Each column as a line
setLineWidth(vThetaHd,1.5);
fontName='Times New Roman';fontWeight='bold';
title(['Frequency=',f0,'MHz'],'FontName',fontName,'fontWeight',fontWeight);
xlabel('Theta(Deg)','FontName',fontName,'fontWeight',fontWeight);
ylabel('RCS(dBsm)','FontName',fontName,'fontWeight',fontWeight);
vThetaAx=gca;set(vThetaAx,'FontName',fontName,'fontWeight',fontWeight);
legend(vThetaComments(1:2:end-1),'location','bestoutside');
%     vThetaAx.FontName='Times New Roman';vThetaAx.FontWeight='bold';
subplot(2,2,2);
vThetaHd=plot(phiVec,rcsDB');%Each column as a line
setLineWidth(vThetaHd,1.5);
fontName='Times New Roman';fontWeight='bold';
title(['Frequency=',f0,'MHz'],'FontName',fontName,'fontWeight',fontWeight);
xlabel('Phi(Deg)','FontName',fontName,'fontWeight',fontWeight);
ylabel('RCS(dBsm)','FontName',fontName,'fontWeight',fontWeight);
vThetaAx=gca;set(vThetaAx,'FontName',fontName,'fontWeight',fontWeight);
legend(vPhiComments(1:2:end-1),'location','bestoutside');
%     vThetaAx.FontName='Times New Roman';vThetaAx.FontWeight='bold';
subplot(2,2,3);
[thetaGrid,phiGrid]=ndgrid(thetaVec,phiVec);
[picHd,axHd,cbHd]=plot3D(thetaGrid,phiGrid,rcs,rcs);
title(['Frequency=',f0,'MHz'],'fontName',fontName,'fontWeight',fontWeight);
xlabel('Theta(deg)','fontName',fontName,'fontWeight',fontWeight);
ylabel('Phi(deg)','fontName',fontName,'fontWeight',fontWeight);
set(axHd,'fontName',fontName,'fontWeight',fontWeight);
cbHd.Label.String='m^2';set(cbHd,'fontName',fontName,'fontWeight',fontWeight);

subplot(2,2,4);
[picHdDB,axHdDB,cbHdDB]=plot3D(thetaGrid,phiGrid,rcsDB,rcsDB);
title(['Frequency=',f0,'MHz'],'fontName',fontName,'fontWeight',fontWeight);
xlabel('Theta(deg)','fontName',fontName,'fontWeight',fontWeight);
ylabel('Phi(deg)','fontName',fontName,'fontWeight',fontWeight);
set(axHdDB,'fontName',fontName,'fontWeight',fontWeight);
cbHdDB.Label.String='dBsm';set(cbHdDB,'fontName',fontName,'fontWeight',fontWeight);

%     gcf.OuterPosition=get(0,'screensize');
set(gcf,'outerposition',get(0,'screensize'));
%     figHd=gcf;
frame=getframe(figHd);
imwrite(frame.cdata,[ExportDir,simName,'(2D_3D).bmp']);

%     count=count+1;
fprintf('\n## %d of %d Done.\n',count,nTxt);
% end
save([ffeDir,matName],'varsMat','varsListOb','varsVec','varsVecOb','rcs_ND','rcsDB_ND','thetaVec','phiVec','-ascii');
disp('Program Done!');
sound(sin(2*pi*25*(1:4000)/500));
