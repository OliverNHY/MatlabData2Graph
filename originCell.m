function dataCell=originCell(varNames,varUnits,varComments,xVec,data)
% originCell({'Theta','RCS'},{'Deg','dB'},{'LongNames','Phi=x'},xVec,data);
% Return a cell with Origin Format.
[~,cs]=size(data);xVec=reshape(xVec,[length(xVec),1]);
% dataCell=cell(rs+3,cs);
indVar=varNames{1};
indVarUnits=varUnits{1};
indVarComments=varComments{1};
nNames=length(varNames)-1;mNameRep=ceil(cs./nNames);
nUnits=length(varUnits)-1;mUnitRep=ceil(cs./nUnits);
nComments=length(varComments)-1;mCommentRep=ceil(cs./nComments);
dpVars=repmat(varNames(2:end),[1,mNameRep]);
dpUnits=repmat(varUnits(2:end),[1,mUnitRep]);
dpComments=repmat(varComments(2:end),[1,mCommentRep]);
x_data_Cell=num2cell([xVec,data]);
dataCell=[indVar,dpVars(1:cs);indVarUnits,dpUnits(1:cs);indVarComments,dpComments(1:cs);x_data_Cell];
end