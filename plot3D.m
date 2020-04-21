function   [picHandle,axHandle,cbHandle]=plot3D(X_grid,Y_grid,Z_grid,color_grid)
%plot3D()  ‰»Îµƒ∂˛Œ¨æÿ’Ûª≠3DÕº

X_max=max(max(X_grid));
Y_max=max(max(Y_grid));
Z_max=max(max(Z_grid));
X_min=min(min(X_grid));
Y_min=min(min(Y_grid));
Z_min=min(min(Z_grid));
if Z_max==Z_min
    Z_max=Z_max+1;
    Z_min=Z_min-1;
end
if Y_max==Y_min
    Y_max=Y_max+1;
    Y_min=Y_min-1;
end
if X_max==X_min
    X_max=X_max+1;
    X_min=X_min-1;
end
surf(X_grid,Y_grid,Z_grid,color_grid);
axis([X_min X_max Y_min Y_max Z_min Z_max]);
shading interp;
colormap(jet);
alpha(0.5);
cbHandle=colorbar;
view([0 0 1]);%z
% view([1 1 0.5]);%3D
% view([0 1 0]);%y
% view([1 0 0]);%x
% view([1 1 1]);
% title(sprintf(' ' ));
fprintf('\n\tPlot3D Done!!!\n');
% axis square
% axis equal
picHandle=gcf;
axHandle=gca;
end

