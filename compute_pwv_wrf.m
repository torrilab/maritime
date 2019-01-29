load('wrf_domain.mat')


Nt     = 1440;
Nz     = 44; 
ggr    = 9.81;

Rgas   = 287; 
cp     = 1004;

p0     = 100000;

root   = '/n/regal/kuang_lab/gtorri/WRFOUT/';
fn     = dir([root,'wrfout_d02*']);
Nf     = length(fn);


rho    = zeros(Nlat,Nlong,Nz);
pwv    = zeros(Nlat,Nlong,Nt);
zi     = zeros(Nlat,Nlong,Nz+1);


for nf = 1
    
    ncp    = netcdf.open([root,fn(nf).name],'NOWRITE');
    
%     ncp    = netcdf.open('/n/regal/kuang_lab/gtorri/WRFOUT/wrfout_d02_2012-03-01_00:00:00', 'NOWRITE');
    
    tic
    for t = 1:Nt
    
        varid  = netcdf.inqVarID(ncp,'QVAPOR');
        qv     = netcdf.getVar(ncp, varid, [0 0 0 (t-1)], [Nlat Nlong Nz 1]);
        
        varid  = netcdf.inqVarID(ncp,'T');
        th     = netcdf.getVar(ncp, varid, [0 0 0 (t-1)], [Nlat Nlong Nz 1])+300;
        
        varid  = netcdf.inqVarID(ncp,'PB');
        p      = netcdf.getVar(ncp, varid, [0 0 0 (t-1)], [Nlat Nlong Nz 1]);
        
        varid  = netcdf.inqVarID(ncp,'P');
        p      = netcdf.getVar(ncp, varid, [0 0 0 (t-1)], [Nlat Nlong Nz 1])+p;
        
        varid  = netcdf.inqVarID(ncp,'PHB');
        z      = netcdf.getVar(ncp, varid, [0 0 0 (t-1)], [Nlat Nlong Nz+1 1])/ggr;
        
        varid  = netcdf.inqVarID(ncp,'PH');
        z      = netcdf.getVar(ncp, varid, [0 0 0 (t-1)], [Nlat Nlong Nz+1 1])/ggr+z;
        
        dz     = z(:,:,2:Nz+1)-z(:,:,1:Nz);
                
        T      = th.*(p./p0).^(Rgas/cp);
        rho    = p./(Rgas.*T);
        
        pwv(:,:,t) = squeeze( sum(rho.*qv.*dz,3) );
        
        disp([t nf])
        
    end
    toc
    
    netcdf.close(ncp);
    
end

save(['pwv_wrf_',num2str(nf),'_new.mat'])



