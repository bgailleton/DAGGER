% PAth to the folder containing the library built with builder.m
addpath("../DAGGER");
system('ldd ../DAGGER/DAGGERInterface.so');
% topo = clib.DAGGER.daggerFD_double_Int_;

% 
ny = 512;
nx = 512;
dx = 50.;
dy = 50.;
xmin = 0.;
ymin = 0.;

topo = clib.DAGGER.quick_fluvial_topo(5, "periodic_EW");


% daggerFD = clib.DAGGER.daggerFD_double_Int_();
% daggerFD.init(nx,ny,dx,dy,xmin,ymin, "periodic_EW");
% 
% f = figure;
% imshow(reshape(topo, [ny,nx]));
% 
% 
% for i = 1:2
%     ttopo = clibConvertArray(clib.matlab.data.TypedArray_double, topo);
%     filledtopo = daggerFD.compute(topo, true);
%     S = (filledtopo(daggerFD.ix) - filledtopo(daggerFD.ixc))/daggerFD.distances;
%     A = daggerFD.get_DA();
% 
% end
