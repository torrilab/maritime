clear all

load('../data/wrf_hgt.mat')
load('../data/wrf_domain.mat')

set(0, 'defaultaxesfontsize', 19)
set(0, 'defaultlinelinewidth', 2)
load('../data/wrf_domain.mat')
load('../data/distance_from_sumatra.mat')

min_lat  = min(lat);
max_lat  = max(lat);

min_long = min(long);
max_long = max(long);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute elevation gradients %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_hgt_x = zeros(size(hgt));
d_hgt_y = zeros(size(hgt));

for j = 2:size(hgt,2)-1
for i = 2:size(hgt,1)-1

    d_hgt_x(i,j) = (hgt(i+1,j)-hgt(i-1,j))/(2*3);
    d_hgt_y(i,j) = (hgt(i,j+1)-hgt(i,j-1))/(2*3);
    
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute climatology for SuGAr %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('../data/gps_data_2008-2013_alph.mat')

Nf = size(gps_data,1);

coor_sugar   = zeros(Nf,3);
dcycle_sugar = zeros(Nf,8);
num_sugar    = zeros(Nf,8);

for nf = 1:Nf
    
    coor = cell2mat(gps_data(nf,2));
    data = cell2mat(gps_data(nf,3));
    
    Nd   = size(data,1);
    hr   = floor(data(:,5)/3)+1;
    
    for nd = 1:Nd
        dcycle_sugar(nf,hr(nd)) = dcycle_sugar(nf,hr(nd)) + data(nd,8);
        num_sugar(nf,hr(nd))    = num_sugar(nf,hr(nd)) + 1;
    end
    
    coor_sugar(nf,1) = coor(1);
    coor_sugar(nf,2) = coor(2);
    coor_sugar(nf,3) = coor(3);
    
end

dcycle_sugar(num_sugar < 10) = NaN;
dcycle_sugar  = dcycle_sugar./num_sugar;

pwv_avg_sugar = nanmean(dcycle_sugar,2);
pwv_amp_sugar = nanmax(dcycle_sugar,[],2)-nanmin(dcycle_sugar,[],2);

para_sugar    = zeros(3,Nf);
x_sugar       = 2*pi*(0:3:21)/24;
for nf = 1:Nf
    
    % Read the data in a format to feed the regression
    y_sugar = dcycle_sugar(nf,:);
    
    idx  = find(isnan(y_sugar));
    if(~isempty(idx))
        idx_a  = idx+1;
        idx_b  = idx-1;
        idx_a(idx_a > 8) = 1;
        idx_b(idx_b < 1) = 8;
        y_sugar(idx)     = 0.5*(y_sugar(idx_a)+y_sugar(idx_b));
    end
    
    % Do the regression
    yr   = (max(y_sugar)-min(y_sugar));
    ym   = nanmean(y_sugar);
    fit  = @(b,x_sugar)  b(1).*(cos(x_sugar - b(2))) + b(3);
    fcn  = @(b) sum((fit(b,x_sugar) - y_sugar).^2);
    para_sugar(:,nf) = fminsearch(fcn, [yr;  pi;  ym]);
    
end
para_sugar(2,:) = para_sugar(2,:)/(2*pi)*24 + 0*7;


% Impose that amplitude be positive
vec1 = para_sugar(1,:);
vec2 = para_sugar(2,:);

vec2(vec1 < 0) = vec2(vec1 < 0)+12;
vec1 = abs(vec1);
para_sugar(1,:) = vec1;

vec2(vec2 < 0)  = vec2(vec2 < 0)  + 24;
vec2(vec2 > 24) = vec2(vec2 > 24) - 24;
para_sugar(2,:) = vec2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute climatology for NVAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('../data/nvaps_data_1988_2009.mat');
dcycle_nvap = cwv_avg;

lat0       = -89.75;
long0      = -179.75;

lat_nvap   = lat0:0.5:-lat0;
long_nvap  = long0:0.5:-long0;

Nlat_nvap  = length(lat_nvap);
Nlong_nvap = length(long_nvap);

[~,min_j]  = min(abs(lat_nvap - min_lat));
[~,max_j]  = min(abs(lat_nvap - max_lat));
[~,min_i]  = min(abs(long_nvap - min_long));
[~,max_i]  = min(abs(long_nvap - max_long));

para_nvap  = zeros(Nlong_nvap,Nlat_nvap,3);
x_nvap     = 2*pi*(0:6:21)/24;

for j = min_j-10:max_j+10     %1:Nlat_nvap
    for i = min_i-10:max_i+10 %1:Nlong_nvap
        
        % Read the data in a format to feed the regression
        y_nvap = squeeze( dcycle_nvap(i,j,:) )';
        
        idx  = find(isnan(y_nvap));
        if(~isempty(idx))
            idx_a = idx+1;
            idx_b = idx-1;
            idx_a(idx_a > 4) = 1;
            idx_b(idx_b < 1) = 4;
            y_nvap(idx)      = 0.5*(y_nvap(idx_a)+y_nvap(idx_b));
        end
        
        % Do the regression
        yr   = (max(y_nvap)-min(y_nvap));
        ym   = nanmean(y_nvap);
        fit  = @(b,x_nvap)  b(1).*(cos(x_nvap - b(2))) + b(3);
        fcn  = @(b) sum((fit(b,x_nvap) - y_nvap).^2);
        para_nvap(i,j,:) = fminsearch(fcn, [yr;  pi;  ym]);
                
    end
end
para_nvap(:,:,2) = para_nvap(:,:,2)/(2*pi)*24 + 0*3;


% Impose that amplitude be positive
vec1             = para_nvap(:,:,1);
vec2             = para_nvap(:,:,2);

vec2(vec1 < 0)   = vec2(vec1 < 0) + 12;
vec1             = abs(vec1);
para_nvap(:,:,1) = vec1;

vec2(vec2 < 0)   = vec2(vec2 < 0)  + 24;
vec2(vec2 > 24)  = vec2(vec2 > 24) - 24;
para_nvap(:,:,2) = vec2;



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparing transect %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% MREK - SDKL - RNDG - PBLI - PBKR - BABI
%  NaN -  37  -  35  -  27  -  26  -  NaN

dcycle_nvap_comp  = zeros(4,4);
dcycle_sugar_comp = zeros(8,4);

n = 0;
for nf = [37,35,27,26]
    n = n+1;
    [~,j] = min(abs(lat_nvap  - coor_sugar(nf,1)));
    [~,i] = min(abs(long_nvap - coor_sugar(nf,2)));

    dcycle_nvap_comp(:,n)  = dcycle_nvap(i,j,:);
    dcycle_sugar_comp(:,n) = dcycle_sugar(nf,:);
end

save('dcycle_nvap_comp.mat', 'dcycle_nvap_comp', 'dcycle_sugar_comp')

figure('Position', [1 1 800 800])
hold on
plot((0:3:24),cat(1,dcycle_sugar_comp,dcycle_sugar_comp(1,:)))
scatter(0:6:24, [dcycle_nvap_comp(:,1)' dcycle_nvap_comp(1,1)],200,'fill', 'markerfacecolor', [0    0.4470    0.7410])
scatter(0:6:24, [dcycle_nvap_comp(:,2)' dcycle_nvap_comp(1,2)],200,'fill', 'markerfacecolor', [0.8500    0.3250    0.0980])
scatter(0:6:24, [dcycle_nvap_comp(:,3)' dcycle_nvap_comp(1,3)],200,'fill', 'markerfacecolor', [0.9290    0.6940    0.1250])
scatter(0:6:24, [dcycle_nvap_comp(:,4)' dcycle_nvap_comp(1,4)],200,'fill', 'markerfacecolor', [0.4940    0.1840    0.5560])
hold off
grid on
box on
xlim([0 24])
set(gca, 'xtick', 0:3:24)
xlabel('UTC (hrs)')
ylabel('PWV (mm)')
legend('SDKL','RNDG','PBLI','PBKR', 'location', 'northwest')


%%%%%%%%%%%%%%%%%%%%
%%% Print output %%%
%%%%%%%%%%%%%%%%%%%%


figure('Position', [1 1 2000 600])

type = 3;
sb1  = subplot(1,3,1);
hold on
contourf(long_nvap,lat_nvap,squeeze( para_nvap(:,:,type) )',200, 'linestyle', 'none')
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
hold off
box on
xlim([94.5 106.5])
ylim([-6 6])
xlabel('Longitude')
ylabel('Latitude')
cb = colorbar;
caxis([38 58])
title(cb, '(mm)')
sb1.Position = [0.04 0.11 0.240 0.8150];

type = 1;
sb2  = subplot(1,3,2);
hold on
contourf(long_nvap,lat_nvap,squeeze( para_nvap(:,:,type) )',200, 'linestyle', 'none')
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
hold off
box on
xlim([94.5 106.5])
ylim([-6 6])
xlabel('Longitude')
ylabel('Latitude')
cb = colorbar;
caxis([0 3])
title(cb, '(mm)')
sb2.Position = [0.37 0.11 0.240 0.8150];

sb3 = subplot(1,3,3);
type = 2;
hold on
contourf(long_nvap,lat_nvap,squeeze( para_nvap(:,:,type) )',200, 'linestyle', 'none')
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
hold off
box on
xlim([94.5 106.5])
ylim([-6 6])
xlabel('Longitude')
ylabel('Latitude')
cb = colorbar;
caxis([0 24])
title(cb, '(UTC)')
set(cb, 'xTick', 0:3:24)
sb3.Position = [0.7 0.11 0.240 0.8150];





vec = zeros(1,Nf);
for nf = 1:Nf
    [~,j] = min(abs(lat  - coor_sugar(nf,1)));
    [~,i] = min(abs(long - coor_sugar(nf,2)));
%     vec(nf) = max( abs(d_hgt_x(i,j)),abs(d_hgt_y(i,j)) );
    vec(nf)  = dist_mat(i,j);
%     vec(nf) = coor_sugar(nf,3);
%     vec(nf) = coor_sugar(nf,1);
end



figure('Position', [1 1 2000 1200])

type = 3;
sb1 = subplot(2,3,1);
hold on
contourf(long_nvap,lat_nvap,squeeze( para_nvap(:,:,type) )',200, 'linestyle', 'none')
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
hold off
box on
xlim([94.5 106.5])
ylim([-6 6])
caxis([38 58])
xlabel('Longitude')
ylabel('Latitude')
cb2 = colorbar;
title(cb2, '(mm)')
sb1.Position = [0.04 0.56 0.240 0.410];

type = 1;
sb2 = subplot(2,3,2);
hold on
contourf(long_nvap,lat_nvap,squeeze( para_nvap(:,:,type) )',200, 'linestyle', 'none')
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
hold off
box on
xlim([94.5 106.5])
ylim([-6 6])
caxis([0 3])
xlabel('Longitude')
ylabel('Latitude')
cb2 = colorbar;
title(cb2, '(mm)')
sb2.Position = [0.37 0.56 0.240 0.410];

type = 2;
sb3 = subplot(2,3,3);
hold on
contourf(long_nvap,lat_nvap,squeeze( para_nvap(:,:,type) )',200, 'linestyle', 'none')
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
for nf = 1:Nf
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
hold off
box on
xlim([94.5 106.5])
ylim([-6 6])
xlabel('Longitude')
ylabel('Latitude')
caxis([0 24])
cb1 = colorbar;
set(cb1,'YTick',0:3:24)
title(cb1, '(UTC)')
sb3.Position = [0.7 0.56 0.240 0.410];


sb4 = subplot(2,3,4);
type = 3;
hold on
line([0 100],[0 100], 'color', 'k')
for nf = 1:Nf
    [~,i] = min(abs(long_nvap-coor_sugar(nf,2)));
    [~,j] = min(abs(lat_nvap-coor_sugar(nf,1)));
    if(~isnan(para_sugar(3,nf)))
    scatter(para_sugar(type,nf), para_nvap(i,j,type), 150, vec(nf), 'fill') 
    end
end
hold off
grid on
box  on
% xlim([46 56])
% ylim([46 56])
xlim([38 58])
ylim([38 58])
% set(gca,'xtick',46:2:58)
% set(gca,'ytick',46:2:58)
set(gca,'xtick',38:4:58)
set(gca,'ytick',38:4:58)
xlabel('SuGAr mean (mm)')
ylabel('NVAP-M mean (mm)')
cb4 = colorbar;
title(cb4, '(km)')
caxis([0 100])
caxis([-50 150])
sb4.Position = [0.04 0.07 0.240 0.410];

sb5 = subplot(2,3,5);
type = 1;
hold on
line([0 100],[0 100], 'color', 'k')
for nf = 1:Nf
    [~,i] = min(abs(long_nvap-coor_sugar(nf,2)));
    [~,j] = min(abs(lat_nvap-coor_sugar(nf,1)));
    if(~isnan(para_sugar(3,nf)))
    scatter(para_sugar(type,nf), para_nvap(i,j,type), 150, vec(nf), 'fill') 
    end
end
hold off
grid on
box  on
xlim([0 3])
ylim([0 3])
xlabel('SuGAr amplitude (mm)')
ylabel('NVAP-M amplitude (mm)')
cb5 = colorbar;
title(cb5, '(km)')
caxis([0 100])
caxis([-50 150])
sb5.Position = [0.37 0.07 0.240 0.410];

sb6 = subplot(2,3,6);
type = 2;
hold on
line([0 100],[0 100], 'color', 'k')
for nf = 1:Nf
    [~,i] = min(abs(long_nvap-coor_sugar(nf,2)));
    [~,j] = min(abs(lat_nvap-coor_sugar(nf,1)));
    if(~isnan(para_sugar(3,nf)))
    scatter(para_sugar(type,nf), para_nvap(i,j,type), 150, vec(nf), 'fill') 
    end
end
hold off
grid on
box  on
set(gca,'xtick',0:3:24)
set(gca,'ytick',0:3:24)
xlim([0 24])
ylim([0 24])
xlabel('SuGAr phase (UTC)')
ylabel('NVAP-M phase (UTC)')
cb6 = colorbar;
title(cb6, '(km)')
caxis([0 100])
caxis([-50 150])
sb6.Position = [0.7 0.07 0.240 0.410];



