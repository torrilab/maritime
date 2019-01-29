%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Accumulate data and compute averages  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'DefaultAxesFontSize', 19)
clear all
load('../data/gps_data_2008-2013_alph.mat')

load('../data/wrf_hgt.mat')
load('../data/wrf_domain.mat')
load('../data/distance_from_sumatra.mat')


Nf         = size(gps_data,1);
coor_gps   = zeros(Nf,3);
pwv_avg_dc = zeros(Nf,8);
num_avg_dc = zeros(Nf,8);
dist_gps   = zeros(Nf,1);

for nf = 1:Nf
    
    coor = cell2mat(gps_data(nf,2));
    coor_gps(nf,1) = coor(1);
    coor_gps(nf,2) = coor(2);
    coor_gps(nf,3) = coor(3);
    
    [~,i] = min(abs(long-coor_gps(nf,2)));
    [~,j] = min(abs(lat-coor_gps(nf,1)));
    
    dist_gps(nf) = dist_mat(i,j);
    
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





% Correlate height with mean PWV

figure('Position', [1 1 1300 650])

subplot(1,2,1)
scatter(para(3,:), coor_gps(:,3)', 200, 'fill')

ylabel('Altitude (m)')
xlabel('Mean - SuGAr (mm)')
ylim([-20 1000])
xlim([39 55])
set(gca, 'xtick', 39:4:55)
set(gca, 'ytick', 0:200:1000)
grid on
box on


subplot(1,2,2)
scatter(dist_gps(:)',para(1,:), 200, 'fill')

ylabel('Amplitude - SuGAr (mm)')
xlabel('Distance from coast (km)')
set(gca,'xdir','reverse')
xlim([-50 150])
grid on
box on




scatter(dist_gps(:)',para(2,:), 200, 'fill')

ylabel('Amplitude - SuGAr (mm)')
xlabel('Distance from coast (km)')
set(gca,'xdir','reverse')
xlim([-50 150])
grid on
box on

