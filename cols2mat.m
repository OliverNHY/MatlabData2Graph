% cols=[repmat((1:10)',10,1),reshape(repmat((11:20),10,1),100,1),(1:100)'];
% [zMat,vecs]=cols2mat1(cols)
function [zMat,vecs]=cols2mat(cols)
% cols is column matrix, in which the last columns is depended variate
% size(zMat)=[length(zMat(:,1)),length(zMat(:,2)),...,,length(zMat(:,length(cols(1,:))-1))]
% Support cols is 3 columns in current version
[nr,nc]=size(cols);
vecs=cell(1,nc-1);
for icol=1:nc-1
    vecs{icol}=unique(cols(:,icol));
    vecL(icol)=length(vecs{icol});
end
if nc~=3
    fprintf('\ncols2mat(): Should be 3 columns.\n');
else
    zMat=NaN(vecL(1),vecL(2));
    for index=1:nr
        z_ir=vecs{1}==cols(index,1);
        z_ic=vecs{2}==cols(index,2);
        zMat(z_ir,z_ic)=cols(index,end);
    end
end
end

