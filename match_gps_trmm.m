


% Load the GPS stations
gps_name = 'SDKL_2012-2012';
gps_coo  = [2.796642489    98.192918679     894.4091617];

root     = '/n/kuanglfs/gtorri/';
load([root,'data_GPS/',gps_name,'_postproc.mat'])


Np = size(gps,1);

% Read the TRMM data

root = '/n/regal/kuang_lab/gtorri/data_TRMM/3-hourly/';

pr_v = zeros(1,Np);
flag = 0;

for np = 1:Np
    
    yr = gps(np,1);
    mt = gps(np,3);
    dy = gps(np,4);
    hr = gps(np,5);
    
    filename = ['3B42.',num2str(yr),num2str(mt,'%02d'),num2str(dy,'%02d'),'.',num2str(hr,'%02d'),'.7.SUB.nc'];
    
    f = dir([root,filename]);
    if(size(f,1) > 0)
        
        % Open file if it exists
        nc = netcdf.open([root,filename], 'NOWRITE');
        
        if(flag == 0)
            
            % The first time, read the coordinates
            varid = netcdf.inqVarID(nc,'latitude');
            lat   = netcdf.getVar(nc, varid);
            
            varid = netcdf.inqVarID(nc,'longitude');
            lon   = netcdf.getVar(nc, varid);
            
            Nlat  = length(lat);
            Nlon  = length(lon);
            
            [~,jy] = (min( abs(lat-gps_coo(1)) ));
            [~,ix] = (min( abs(lon-gps_coo(2)) ));
            
            flag = 1;
            
        end
        
        varid = netcdf.inqVarID(nc,'pcp');
        prec  = netcdf.getVar(nc, varid);
        
        pr_v(np) = prec(ix,jy);
        
    end

    
end

gps  = cat(2,gps,pr_v');

root = '/n/kuanglfs/gtorri/';
save([root,'data_GPS/',gps_name,'_gps.mat'], 'gps')
