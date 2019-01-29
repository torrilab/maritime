root    = '/n/kuanglfs/gtorri/NVAPS/';


cwv_avg = zeros(720,360,4);
cwv_std = zeros(720,360,4);
cwv_num = zeros(720,360,4);

yr_i    = 1988;
yr_f    = 2009;


for year = yr_i:yr_f
    for day = 1:366
        
        file = [root,num2str(year),'.',num2str(day,'%03.f'),'.NVAPS.tpw.4xdaily.nc'];
        
        if exist(file,'file') == 2
            
            nc = netcdf.open(file, 'NOWRITE');
           
            varid   = netcdf.inqVarID(nc,'water_vapor');
            cwv_tmp = squeeze(netcdf.getVar(nc, varid));
            
            id      = find(cwv_tmp > -9);
            [I,J,T] = ind2sub([720,360,4],id);
            
            for n = 1:length(I)
                cwv_avg(I(n),J(n),T(n)) = cwv_avg(I(n),J(n),T(n)) + cwv_tmp(I(n),J(n),T(n));
                cwv_num(I(n),J(n),T(n)) = cwv_num(I(n),J(n),T(n)) + 1;
            end
            netcdf.close(nc)
            disp([day year])
            
        end
        
    end
    
end
cwv_avg = cwv_avg./cwv_num;


for year = yr_i:yr_f
    for day = 1:366
        
        file = [root,num2str(year),'.',num2str(day,'%03.f'),'.NVAPS.tpw.4xdaily.nc'];
        
        if exist(file,'file') == 2
            
            nc = netcdf.open(file, 'NOWRITE');
            
            varid   = netcdf.inqVarID(nc,'water_vapor');
            cwv_tmp = squeeze(netcdf.getVar(nc, varid));
            
            id      = find(cwv_tmp > -9 & cwv_avg > 0);
            [I,J,T] = ind2sub([720,360,4],id);
            
            for n = 1:length(I)
                if(cwv_avg(I(n),J(n),T(n)) > 0)
                    cwv_std(I(n),J(n),T(n)) = cwv_std(I(n),J(n),T(n)) + ...
                        (cwv_tmp(I(n),J(n),T(n)) - cwv_avg(I(n),J(n),T(n)))^2 ;
                end
            end
            netcdf.close(nc)
            disp([day year])
            
        end
        
    end
    
end
cwv_std = sqrt( cwv_std./(cwv_num-1) );


save(['nvaps_data_',num2str(yr_i),'_',num2str(yr_f),'.mat'])
