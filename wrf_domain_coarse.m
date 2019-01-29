

Nlat   = 444;
Nlong  = 777;

root   = '/n/regal/kuang_lab/gtorri/';
fn     = dir([root,'wrfout_d01*']);
Nf     = length(fn);

nf     = 1;

ncp    = netcdf.open(fn(nf).name,'NOWRITE');


varid  = netcdf.inqVarID(ncp,'XLAT');
lat    = netcdf.getVar(ncp, varid, [0 0 0], [1 Nlat 1]);

varid  = netcdf.inqVarID(ncp,'XLONG');
long   = netcdf.getVar(ncp, varid, [0 0 0], [Nlong 1 1]);

varid  = netcdf.inqVarID(ncp,'HGT');
hgt    = netcdf.getVar(ncp, varid, [0 0 0], [Nlong Nlat 1]);

varid  = netcdf.inqVarID(ncp,'LANDMASK');
lmask  = netcdf.getVar(ncp, varid, [0 0 0], [Nlong Nlat 1]);


netcdf.close(ncp);


fac     = 4;
nlong_c = floor(Nlong/fac);
nlat_c  = floor(Nlat/fac);

lat_c   = zeros(1,nlat_c);
long_c  = zeros(1,nlong_c);
lmask_c = zeros(nlong_c,nlat_c);



for j = 1:nlat_c
    for n2 = 1:fac
        lat_c(j) = lat_c(j) + lat(fac*(j-1)+n2)/fac;
    end
end



for i = 1:nlong_c
    for n1 = 1:fac
        long_c(i) = long_c(i) + long(fac*(i-1)+n1)/fac;
    end
end



for j = 1:nlat_c
for i = 1:nlong_c
        
    for n1 = 1:fac
    for n2 = 1:fac
       lmask_c(i,j) = lmask_c(i,j) + lmask(fac*(i-1)+n1,fac*(j-1)+n2)/fac^2;
    end
    end
        
end
end



save('wrf_domain_coarse.mat', 'nlat_c', 'nlong_c', 'lat_c',  'long_c', 'lmask_c')



