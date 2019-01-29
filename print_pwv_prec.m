


% Load the GPS stations
gps_name = 'LEWK_2010-2013';

load([gps_name,'_gps.mat'])


gps_coo  = [2.923604066    95.804053163      12.9097333];


lat_in = gps_coo(1);
lon_in = gps_coo(2);
d_meth = 'fast';
d_max  = 1000*1000;
dist   = (dist_from_coast(lat_in,lon_in,d_meth,d_max))/1000; % in km

LAND   = landmask(gps_coo(1),gps_coo(2));



max_pwv = 70;
d_pwv   = 2;
n_pwv   = floor(max_pwv/d_pwv);

fit_v = zeros(1,n_pwv);
num_v = zeros(1,n_pwv);

for nh = 1:n_pwv
    
    id = find(gps(:,8) > d_pwv*(nh-1) & gps(:,8) <= d_pwv*nh);
    
    for nid = 1:length(id)
        
        fit_v(nh) = fit_v(nh) + gps(id(nid),9);
        num_v(nh) = num_v(nh) + 1;
        
    end
    
end







