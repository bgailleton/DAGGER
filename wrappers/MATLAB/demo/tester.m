addpath("../DAGGER");

ny = 500;
nx = 500;
dx = 50.;
dy = 50.;
xmin = 0.;
ymin = 0.;

topo = rand(nx * ny, 1);

daggerFD = clib.DAGGER.daggerFD_double_Int_();
daggerFD.init(nx,ny,dx,dy,xmin,ymin, "periodic_EW");

f = figure;
imshow(reshape(topo, [ny,nx]));


for i = 1:2
    ttopo = clibConvertArray(clib.matlab.data.TypedArray_double, topo);
    filledtopo = daggerFD.compute(topo, true);
    S = (filledtopo(daggerFD.ix) - filledtopo(daggerFD.ixc))/daggerFD.distances;
    A = daggerFD.get_DA();

end
