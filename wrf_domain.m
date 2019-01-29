

Nlat   = 402;
Nlong  = 402;

root   = '/n/regal/kuang_lab/gtorri/';
fn     = dir([root,'wrfout_d02*']);
Nf     = length(fn);

nf     = 1;

ncp    = netcdf.open(fn(nf).name,'NOWRITE');


varid  = netcdf.inqVarID(ncp,'XLAT');
lat    = netcdf.getVar(ncp, varid, [0 0 0], [1 Nlat 1]);

varid  = netcdf.inqVarID(ncp,'XLAT_V');
lat_i  = netcdf.getVar(ncp, varid, [0 0 0], [1 Nlat+1 1]);


varid  = netcdf.inqVarID(ncp,'XLONG');
long   = netcdf.getVar(ncp, varid, [0 0 0], [Nlong 1 1]);

varid  = netcdf.inqVarID(ncp,'XLONG_U');
long_i = netcdf.getVar(ncp, varid, [0 0 0], [Nlong+1 1 1]);

varid  = netcdf.inqVarID(ncp,'HGT');
hgt    = netcdf.getVar(ncp, varid, [0 0 0], [Nlat Nlong 1]);

varid  = netcdf.inqVarID(ncp,'LANDMASK');
lmask  = netcdf.getVar(ncp, varid, [0 0 0], [Nlat Nlong 1]);


netcdf.close(ncp);

save('wrf_domain.mat', 'Nlat', 'Nlong', 'lat', 'lat_i', 'long', 'long_i', 'lmask', 'hgt')

