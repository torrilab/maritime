clear all
set(0, 'DefaultAxesFontSize', 19)

gps_name = 'LEWK';

yr_i   = 2010;
yr_f   = 2013;
N_data = 0;
gps    = zeros(0,8);

for yr = yr_i:yr_f
    filename = [gps_name,'_',num2str(yr), '_30_hour_GIPSY_PWV_met_data.dat'];
    
    gps_raw = importdata(['../data/gps/',filename]);
    pwv_raw = squeeze(gps_raw(:,9));
    id_r    = find(pwv_raw > 0);
    N_tmp   = length(id_r);
    gps_tmp = zeros(N_tmp,7);
    
    for m = 1:N_tmp
        n = id_r(m);
        
        gps_tmp(m,1) = gps_raw(n,1); % Year
        gps_tmp(m,2) = gps_raw(n,2); % Time of year
        gps_tmp(m,3) = gps_raw(n,3); % Month
        gps_tmp(m,4) = gps_raw(n,4); % Day
        gps_tmp(m,5) = gps_raw(n,5); % Hour
        gps_tmp(m,6) = gps_raw(n,7); % Pressure
        gps_tmp(m,7) = gps_raw(n,8); % Temperature
        gps_tmp(m,8) = gps_raw(n,9); % PWV
        
    end
    
    N_data = N_data + N_tmp;
    gps = cat(1,gps_tmp,gps);
    
end



% Diurnal cycle

dcyc_mean = zeros(24,12); 
dcyc_std  = zeros(24,12); 
norm      = zeros(24,12); 

for m = 1:12
    
    id_m = find(gps(:,3) == m);
    
    for p = 1:length(id_m)
        
        n = id_m(p);
        hr  = gps(n,5)+1;
        pwv = gps(n,8);
        dcyc_mean(hr,m) = dcyc_mean(hr,m) + gps(n,8);
        norm(hr,m) = norm(hr,m)+1;
        
    end
    
    dcyc_mean(:,m) = dcyc_mean(:,m)./norm(:,m);
    
    for p = 1:length(id_m)
    
        n = id_m(p);
        hr  = gps(n,5)+1;
        pwv = gps(n,8);
        dcyc_std(hr,m) = dcyc_std(hr,m) + (gps(n,8)-dcyc_mean(hr,m)).^2;
    
    end
    dcyc_std(:,m) = sqrt(dcyc_std(:,m)./(norm(:,m)-1));
    
    
end
dcyc_stdm = dcyc_std./sqrt(norm);



dcyc_mean_tot = zeros(1,24); 
dcyc_std_tot  = zeros(1,24); 
norm_tot      = zeros(1,24); 

for n = 1:N_data
    
    hr  = gps(n,5)+1;
    dcyc_mean_tot(hr) = dcyc_mean_tot(hr) + gps(n,8);
    norm_tot(hr) = norm_tot(hr) + 1;
    
end

dcyc_mean_tot = dcyc_mean_tot./norm_tot;
dcyc_mean_tot(norm_tot < 100) = NaN;

for n = 1:N_data
    
    hr  = gps(n,5)+1;
    dcyc_std_tot(hr) = dcyc_std_tot(hr) + (gps(n,8) - dcyc_mean_tot(hr)).^2;
    
end
dcyc_std_tot = sqrt(dcyc_std_tot./(norm_tot-1));
dcyc_stdm_tot = dcyc_std_tot./sqrt(norm_tot);

vec_1(1:8) = dcyc_mean_tot(1:3:24);
vec_2(1:8) = dcyc_stdm_tot(1:3:24);

figure('Position', [1 1 800 800])
errorbar((0:7)*3,vec_1,vec_2)
grid on
box on
set(gca,'xTick', 0:3:21)
xlabel('Time (hr)')
ylabel('PWV (mm)')
xlim([0 21])
title([gps_name, ' 2011-2013'])


save([gps_name,'_',num2str(yr_i),'-',num2str(yr_f),'_mean_std.mat'])

