% load('wrf_domain.mat')
% 
% 
% pwv_cosmic = zeros(Nlat,Nlong,8);
% num_cosmic = zeros(Nlat,Nlong,8);
% 
% for year = 2010:2012
%     
%     data_raw = importdata(['../data/cosmic_3hr/',num2str(year),'.txt']);
%     Ndata    = size(data_raw.data,1);
%     
%     lat_v    = data_raw.data(:,5);
%     long_v   = data_raw.data(:,6);
%     
%     id_list  = find(lat_v >= min(lat) & lat_v <= max(lat) & long_v >= min(long) & long_v <= max(long));
%     
%     for m = 1:length(id_list)
%         
%         nd = id_list(m);
%         
%         lat_tmp  = data_raw.data(nd,5);
%         [~,i] = min(abs( lat_i-lat_tmp ));
%         if(lat_i(i) >= lat_tmp)
%             i = i-1;
%         end
%         
%         long_tmp = data_raw.data(nd,6);
%         [~,j] = min(abs( long_i-long_tmp ));
%         if(long_i(j) >= long_tmp)
%             j = j-1;
%         end
%         
%         t   = floor(data_raw.data(nd,3)/3) + 1;
%         pwv = data_raw.data(nd,7);
%         
%         pwv_cosmic(i,j,t) = pwv_cosmic(i,j,t)+pwv;
%         num_cosmic(i,j,t) = num_cosmic(i,j,t)+1;
%         
%     end
%     
%     disp(year)
%     
% end
% 




% Load SuGAr data
load('../data/wrf_domain.mat')
load('../data/gps_data_2008-2013_alph.mat')

N_gps    = size(gps_data,1);
coo_list = cell2mat(gps_data(:,2));

min_dist = 1.5; %0.25;

sugar_cosmic = zeros(0,2);


for year = 2008:2013
    
    data_raw = importdata(['../data/cosmic_3hr/',num2str(year),'.txt']);
    Ndata    = size(data_raw.data,1);
    
    lat_v    = data_raw.data(:,5);
    long_v   = data_raw.data(:,6);
    
    id_list  = find(lat_v >= min(lat) & lat_v <= max(lat) & long_v >= min(long) & long_v <= max(long));
    
    for m = 1:length(id_list)
        
        nd = id_list(m);
        
        lat_tmp  = lat_v(nd);
        long_tmp = long_v(nd);
        dist = sqrt( (lat_tmp-coo_list(:,1)).^2 + (long_tmp-coo_list(:,2)).^2 );
        
        [val,pos] = min(dist);
        
        if(val < min_dist)
            
            dy = data_raw.data(nd,2);
            hr = floor(data_raw.data(nd,3)/3);
            
            sugar_data = cell2mat(gps_data(pos,3));
            sugar_yr   = sugar_data(:,1);
            sugar_dy   = datenum(year,sugar_data(:,3),sugar_data(:,4)) - datenum(year,1,1);
            sugar_hr   = sugar_data(:,5);
            
            idx = find(sugar_yr == year & sugar_dy == dy & sugar_hr == hr);
            
            if(~isempty(idx))
                sugar_pwv  = sugar_data(idx,8);
                cosmic_pwv = data_raw.data(nd,7);
                sugar_cosmic = cat(1,sugar_cosmic,[sugar_pwv,cosmic_pwv]);
            end
            
            
        end
        
    end
    
    disp(year)
    
end


X = [ones(size(sugar_cosmic,1),1) sugar_cosmic(:,1)];
y = sugar_cosmic(:,2);

figure('Position', [1 1 650 650])
hold on
line([0 100], [0 100], 'color', [0.8500    0.3250    0.0980], 'linewidth', 2)
scatter(sugar_cosmic(:,1),sugar_cosmic(:,2), 250,'fill', 'markeredgecolor', 'k')
hold off
grid on
box on

xlim([35 65])
ylim([35 65])

xlabel('PWV from SuGAr (mm)')
ylabel('PWV from COSMIC (mm)')




