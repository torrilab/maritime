%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Accumulate data and compute averages  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'DefaultAxesFontSize', 19)
clear all
load('../data/gps_data_2008-2013_alph.mat')

load('../data/wrf_hgt.mat')
load('../data/wrf_domain.mat')


Nf         = size(gps_data,1);
coor_gps   = zeros(Nf,3);
pwv_avg_dc = zeros(Nf,8);
num_avg_dc = zeros(Nf,8);


for nf = 1:Nf
    
    coor = cell2mat(gps_data(nf,2));
    coor_gps(nf,1) = coor(1);
    coor_gps(nf,2) = coor(2);
    coor_gps(nf,3) = coor(3);
    
    data = cell2mat(gps_data(nf,3));
    Nd   = size(data,1);
    
    hr   = data(:,5);
    for nd = 1:Nd
        tmp = floor(hr(nd)/3)+1;
        pwv_avg_dc(nf,tmp) = pwv_avg_dc(nf,tmp) + data(nd,8);
        num_avg_dc(nf,tmp) = num_avg_dc(nf,tmp) + 1;
    end
    
end
num_avg_dc(num_avg_dc < 10) = NaN;
tmp_1 = nansum(pwv_avg_dc,2);
tmp_2 = nansum(num_avg_dc,2);
vec   = (tmp_1./tmp_2);

pwv_avg_dc = pwv_avg_dc./num_avg_dc;



figure('Position', [1 1 800 800])
hold on
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),100,vec(nf),'fill');
end
hold off
box on
xlim([94 110])
ylim([-8 8])
xlabel('Longitude')
ylabel('Latitude')
cb2 = colorbar;
title(cb2, '(mm)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Regression on the climatology %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


para = zeros(3,Nf);
x    = 2*pi*(0:3:21)/24;

for nf = 1:Nf
    
    % Read the data in a format to feed the regression
    y    = pwv_avg_dc(nf,:);
    
    idx  = find(isnan(y));
    if(~isempty(idx))
        
        idx_a  = idx+1;
        idx_b  = idx-1;
        
        idx_a(idx_a > 8) = 1;
        idx_b(idx_b < 1) = 8;
        
        y(idx) = 0.5*(y(idx_a)+y(idx_b));
        
    end
    
    yr   = (max(y)-min(y));
    
    % Do the regression
    ym   = nanmean(y);
    
    fit  = @(b,x)  b(1).*(cos(x - b(2))) + b(3);
    fcn  = @(b) sum((fit(b,x) - y).^2);
    para(:,nf) = fminsearch(fcn, [yr;  pi;  ym]);
    
end
para(2,:) = para(2,:)/(2*pi)*24 + 0*7;


% Impose that amplitude be positive 
vec1 = para(1,:);
vec2 = para(2,:);

vec2(vec1 < 0) = vec2(vec1 < 0)+12;
vec1 = abs(vec1);
para(1,:) = vec1;

vec2(vec2 < 0)  = vec2(vec2 < 0)  + 24; 
vec2(vec2 > 24) = vec2(vec2 > 24) - 24; 
para(2,:) = vec2;






figure('Position', [1 1 650 650])
hold on
% geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
contour(long,lat,hgt', 500:500:4000, 'linewidth', 2, 'color','k')
contour(long,lat,hgt', 1:100:400, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),250,para(2,nf),'fill', 'markeredgecolor', 'k');
end
hold off
box on
xlim([95 106]) %([94 110])
ylim([-6 6]) %([-8 8])
xlabel('Longitude')
ylabel('Latitude')
caxis([0 24])
cb1 = colorbar;
set(cb1,'YTick',0:3:24)
title(cb1, '(UTC)')
title('\fontsize{35}Phase')


figure('Position', [1 1 650 650])
hold on
% geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
contour(long,lat,hgt', 500:500:4000, 'linewidth', 2, 'color','k')
contour(long,lat,hgt', 1:100:400, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),250,para(1,nf),'fill', 'markeredgecolor', 'k');
end
hold off
box on
xlim([95 106]) %([94 110])
ylim([-6 6]) %([-8 8])
caxis([0 3])
xlabel('Longitude')
ylabel('Latitude')
cb2 = colorbar;
title(cb2, '(mm)')
title('\fontsize{35}Amplitude')


figure('Position', [1 1 650 650])
hold on
% geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
contour(long,lat,hgt', 500:500:4000, 'linewidth', 2, 'color','k')
contour(long,lat,hgt', 1:100:400, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),300,para(3,nf),'fill', 'markeredgecolor', 'k');
end
hold off
box on
xlim([95 106]) %([94 110])
ylim([-6 6]) %([-8 8])
caxis([38 58])
xlabel('Longitude')
ylabel('Latitude')
cb2 = colorbar;
title(cb2, '(mm)')
title('\fontsize{35}Mean')





figure('Position', [1 1 2000 600])

sb1 = subplot(1,3,1);
hold on
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
% contour(long,lat,hgt', 1:100:400, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),300,para(3,nf),'fill', 'markeredgecolor', 'k');
end
hold off
box on
xlim([94 106]) %([94 110])
ylim([-6 6]) %([-8 8])
caxis([38 58])
xlabel('Longitude')
ylabel('Latitude')
cb2 = colorbar;
title(cb2, '(mm)')
% title('\fontsize{35}Mean')
sb1.Position = [0.04 0.11 0.240 0.8150];


sb2 = subplot(1,3,2);
hold on
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
% contour(long,lat,hgt', 1:100:400, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),250,para(1,nf),'fill', 'markeredgecolor', 'k');
end
hold off
box on
xlim([94 106]) %([94 110])
ylim([-6 6]) %([-8 8])
caxis([0 3])
xlabel('Longitude')
ylabel('Latitude')
cb2 = colorbar;
title(cb2, '(mm)')
% title('\fontsize{35}Amplitude')
sb2.Position = [0.37 0.11 0.240 0.8150];


sb3 = subplot(1,3,3);
hold on
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
% contour(long,lat,hgt', 1:100:400, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),250,para(2,nf),'fill', 'markeredgecolor', 'k');
end
hold off
box on
xlim([94 106]) %([94 110])
ylim([-6 6]) %([-8 8])
xlabel('Longitude')
ylabel('Latitude')
caxis([0 24])
cb1 = colorbar;
set(cb1,'YTick',0:3:24)
title(cb1, '(UTC)')
% title('\fontsize{35}Phase')
sb3.Position = [0.7 0.11 0.240 0.8150];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Regression on all the data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


para_all = zeros(3,Nf);

for nf = 1:Nf

    data = cell2mat(gps_data(nf,3));
    
    x    = 2*pi*((data(:,1)-2010).*yeardays(data(:,1)) + data(:,2));
    [xs,I] = sort(x);
    
    % Read the data in a format to feed the regression
    y    = data(:,8);
    y    = y(I);
    yr   = (max(y)-min(y));

    % Do the regression
    ym   = mean(y);
    fit  = @(b,x)  b(1).*(cos(x - b(2))) + b(3);
    fcn  = @(b) sum((fit(b,x) - y).^2);
    para_all(:,nf) = fminsearch(fcn, [yr;  0;  ym]);
    
end


% Normalize time lag
para_all(2,:) = para_all(2,:)/(2*pi)*24;


% Impose that amplitude be positive 
vec1 = para_all(1,:);
vec2 = para_all(2,:);

vec2(vec1 < 0) = vec2(vec1 < 0)+12;
vec1 = abs(vec1);

para_all(1,:) = vec1;
para_all(2,:) = vec2;


figure('Position', [1 1 800 800])
hold on
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
for nf = 1:Nf
    scatter(coor_gps(nf,2),coor_gps(nf,1),100,para_all(1,nf),'fill');
end
hold off
box on
xlim([94 110])
ylim([-8 8])
xlabel('Longitude')
ylabel('Latitude')




