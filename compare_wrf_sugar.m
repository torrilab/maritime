clear all

year0 = 2012;

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultaxesfontsize', 19)
branch = '../data/';

load([branch,'wrf_hgt.mat']);
load([branch,'wrf_domain.mat']);
load([branch,'distance_from_sumatra.mat']);
load([branch,'dcycle_wrf_',num2str(year0),'_UTC.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute climatology for SuGAr %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([branch,'gps_data_2008-2013_alph.mat'])

Nsgr = size(gps_data,1);

coor_sugar   = zeros(Nsgr,2);
dcycle_sugar = zeros(Nsgr,8);
num_sugar    = zeros(Nsgr,8);

for nf = 1:Nsgr
    
    coor = cell2mat(gps_data(nf,2));
    data = cell2mat(gps_data(nf,3));
        
%     ind  = find(data(:,1) == 2011 | data(:,1) == 2012);
    
    ind  = find(data(:,1) == year0);

    day  = floor(data(ind,2));
    hr   = floor(data(ind,5)/3)+1;
    pwv  = data(ind,8);
        
    for m = 1:length(ind)
        dcycle_sugar(nf,hr(m)) = dcycle_sugar(nf,hr(m)) + pwv(m);
           num_sugar(nf,hr(m)) =    num_sugar(nf,hr(m)) + 1;
    end
    
    coor_sugar(nf,1) = coor(1);
    coor_sugar(nf,2) = coor(2);
    
end

dcycle_sugar(num_sugar < 10) = NaN;
dcycle_sugar  = dcycle_sugar./num_sugar;

pwv_avg_sugar = nanmean(dcycle_sugar,2);
pwv_amp_sugar = nanmax(dcycle_sugar,[],2)-nanmin(dcycle_sugar,[],2);


para_sugar    = zeros(3,Nsgr);
x_sugar       = 2*pi*(0:3:21)/24;
for nf = 1:Nsgr
    
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
para_sugar(2,:) = para_sugar(2,:)/(2*pi)*24;


% Impose that amplitude be positive
vec1 = para_sugar(1,:);
vec2 = para_sugar(2,:);

vec2(vec1 < 0) = vec2(vec1 < 0)+12;
vec1 = abs(vec1);
para_sugar(1,:) = vec1;

vec2(vec2 < 0)  = vec2(vec2 < 0)  + 24;
vec2(vec2 > 24) = vec2(vec2 > 24) - 24;
para_sugar(2,:) = vec2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare two stations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(1==0)
% LAIS = Nsgr-6
% NGNG = Nsgr-12

dcycle_point_wrf   = zeros(48,2);
dcycle_point_sugar = zeros(8,2);
vec                = [Nsgr-6, Nsgr-12];

% Compare the climatologies
for m = 1:2

    nf = vec(m);
    
    [~,i] = min(abs(long-coor_sugar(nf,2)));
    [~,j] = min(abs(lat-coor_sugar(nf,1)));
        
    for a1 = -1:1
        for a2 = -1:1
            dcycle_point_wrf(:,m) = dcycle_point_wrf(:,m) + ...
                squeeze( dcycle_wrf(i+a1,j+a2,:)/9 );
        end
    end
    
    dcycle_point_sugar(:,m) = dcycle_sugar(nf,:);
    
    disp(nf)
    
end
    

figure('Position', [1 1 800 800])
hold on
plot(0:0.5:23.5, dcycle_point_wrf)
% plot(0:0.5:23.5, tmp1)
% plot(0:0.5:23.5, tmp2)
for n = 1:8
    scatter(3*(n-1),dcycle_point_sugar(n,1), 'fill', 'markerfacecolor', [0    0.4470    0.7410])
    scatter(3*(n-1),dcycle_point_sugar(n,2), 'fill', 'markerfacecolor', [0.8500    0.3250    0.0980])
end
hold off
grid on
box on
xlim([0 24])
set(gca, 'xtick', 0:3:24)
xlabel('Time UTC (hrs)')
ylabel('PWV (mm)')
legend('LAIS', 'NGNG')



% Compare point by point for two stations

yr0 = 2012;
vec = [Nsgr-6,Nsgr-12];
Nd  = yeardays(yr0)*8;

ifn = zeros(1,2);
jfn = zeros(1,2);

sugar_station = NaN(Nd,2);
for m = 1:2
    
    nf    = vec(m);
    name  = cell2mat(gps_data(nf,1));
    data  = cell2mat(gps_data(nf,3));
    
    yr    = data(:,1);
    idx   = find(yr == 2012);
    ind   = data(idx,2)*8 + 1;
    
    sugar_station(ind,m) = data(idx,8);
    
    [~,ifn(m)] = min(abs(long-coor_sugar(nf,2)));
    [~,jfn(m)] = min(abs(lat-coor_sugar(nf,1)));
    
    disp([nf name])
    
end


wrf_station   = zeros(Nd,2);
for ng = 1:13
        
%     load(['../data/pwv_wrf_',num2str(ng),'_new.mat'])

    load([branch,'/WRFOUT/',num2str(year),'/pwv_wrf_',num2str(ng),'_new.mat'])
    
    ti = 1;
    if(ng == 1)
        ti = 481;
    end
    
    tf = 1440;
    if(ng == 13)
        tf = 277;
    end
    
    nsteps = (tf-ti+1)/6;
    for ns = 1:nsteps
        for m = 1:2
            wrf_station((ng-1)*240+(ti-1)/6+(ns),m) = pwv(ifn(m),jfn(m),ti+(ns-1)*6);
            disp([(ng-1)*240+(ti-1)/6+(ns)  ti+(ns-1)*6])

%             for a1 = -1:1
%                 for a2 = -1:1
%                     wrf_station((ng-1)*240+(ti-1)/6+(ns),m) = ...
%                         wrf_station((ng-1)*240+(ti-1)/6+(ns),m) + ...
%                         pwv(ifn(m)+a1,jfn(m)+a2,ti+(ns-1)*6)/9;
%                 end
%             end

        end
    end
    clear pwv;
    
    disp(ng)
    
end
wrf_station(wrf_station == 0) = NaN;
Nsgr = size(gps_data,1);


figure('Position', [1 1 800 800])
m = 2;
hold on
scatter( sugar_station(:,m),wrf_station(:,m) )
line([0 100],[0 100], 'color', [0.8500    0.3250    0.0980])
hold off
xlim([30 70])
ylim([30 70])
grid on
box on
xlabel('SuGAr PWV (mm)')
ylabel('WRF PWV (mm)')
nm = cell2mat(gps_data(vec(m)));
title(nm)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparing transect %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% MREK - SDKL - RNDG - PBLI - PBKR - BABI
%  NaN -  37  -  35  -  27  -  26  -  NaN

dcycle_wrf_comp   = zeros(48,4);
dcycle_sugar_tmp  = zeros(8,4);
dcycle_sugar_comp = zeros(8,4);

n = 0;
for nf = [37,35,27,26]
    n = n+1;
    [~,j] = min(abs(lat  - coor_sugar(nf,1)));
    [~,i] = min(abs(long - coor_sugar(nf,2)));

    dcycle_wrf_comp(:,n)  = dcycle_wrf(i,j,:);
    dcycle_sugar_tmp(:,n) = dcycle_sugar(nf,:);
end

save('dcycle_wrf_comp.mat', 'dcycle_wrf_comp', 'dcycle_sugar_comp')


figure('Position', [1 1 800 800])
hold on

plot((0:0.5:24),cat(1,dcycle_wrf_comp,dcycle_wrf_comp(1,:)))
scatter(0:3:24, [dcycle_sugar_tmp(:,1)' dcycle_sugar_tmp(1,1)],200,'fill', 'markerfacecolor', [0    0.4470    0.7410])
scatter(0:3:24, [dcycle_sugar_tmp(:,2)' dcycle_sugar_tmp(1,2)],200,'fill', 'markerfacecolor', [0.8500    0.3250    0.0980])
scatter(0:3:24, [dcycle_sugar_tmp(:,3)' dcycle_sugar_tmp(1,3)],200,'fill', 'markerfacecolor', [0.9290    0.6940    0.1250])
scatter(0:3:24, [dcycle_sugar_tmp(:,4)' dcycle_sugar_tmp(1,4)],200,'fill', 'markerfacecolor', [0.4940    0.1840    0.5560])
hold off
grid on
box on
xlim([0 24])
ylim([30 60])
set(gca, 'xtick', 0:3:24)
xlabel('UTC (hrs)')
ylabel('PWV (mm)')
legend('SDKL','RNDG','PBLI','PBKR', 'location', 'northwest')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Print output figures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vec = zeros(1,Nsgr);
for nf = 1:Nsgr
    [~,j] = min(abs(lat  - coor_sugar(nf,1)));
    [~,i] = min(abs(long - coor_sugar(nf,2)));
%     vec(nf) = max( abs(d_hgt_x(i,j)),abs(d_hgt_y(i,j)) );
%     vec(nf)  = dist_mat(i,j);
%     vec(nf) = coor_sugar(nf,3);
    vec(nf) = coor_sugar(nf,1);
end



figure('Position', [1 1 2000 1200])

sb1 = subplot(2,3,1);
hold on
type = 3;
contourf(long,lat,squeeze(para_wrf(:,:,type))',200, 'linestyle', 'none')
for nf = 1:Nsgr
if(~isnan(para_sugar(3,nf)))
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
end
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
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

sb2 = subplot(2,3,2);
hold on
type = 1;
contourf(long,lat,squeeze(para_wrf(:,:,type))',200, 'linestyle', 'none')
for nf = 1:Nsgr
if(~isnan(para_sugar(3,nf)))
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
end
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')
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

sb3 = subplot(2,3,3);
hold on
type = 2;
contourf(long,lat,squeeze(para_wrf(:,:,type))',200, 'linestyle', 'none')
for nf = 1:Nsgr
if(~isnan(para_sugar(3,nf)))
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),250,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
end
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')

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
for nf = 1:Nsgr
    [~,i] = min(abs(long-coor_sugar(nf,2)));
    [~,j] = min(abs(lat-coor_sugar(nf,1)));
    if(~isnan(para_sugar(3,nf)))
    scatter(para_sugar(type,nf), para_wrf(i,j,type), 150, vec(nf), 'fill') 
    end
end
hold off
grid on
box  on
xlim([46 56])
ylim([46 56])
set(gca,'xtick',46:2:58)
set(gca,'ytick',46:2:58)
xlabel('SuGAr mean (mm)')
ylabel('WRF mean (mm)')
sb4.Position = [0.04 0.07 0.240 0.410];

sb5 = subplot(2,3,5);
type = 1;
hold on
line([0 100],[0 100], 'color', 'k')
for nf = 1:Nsgr
    [~,i] = min(abs(long-coor_sugar(nf,2)));
    [~,j] = min(abs(lat-coor_sugar(nf,1)));
    if(~isnan(para_sugar(3,nf)))
    scatter(para_sugar(type,nf), para_wrf(i,j,type), 150,vec(nf), 'fill')   
    end
end
hold off
grid on
box  on
xlim([0 3])
ylim([0 3])
xlabel('SuGAr amplitude (mm)')
ylabel('WRF amplitude (mm)')
sb5.Position = [0.37 0.07 0.240 0.410];

sb6 = subplot(2,3,6);
type = 2;
hold on
line([0 100],[0 100], 'color', 'k')
for nf = 1:Nsgr
    [~,i] = min(abs(long-coor_sugar(nf,2)));
    [~,j] = min(abs(lat-coor_sugar(nf,1)));
    if(~isnan(para_sugar(3,nf)))
    scatter(para_sugar(type,nf), para_wrf(i,j,type), 150,vec(nf), 'fill')  
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
ylabel('WRF phase (UTC)')
sb6.Position = [0.7 0.07 0.240 0.410];














% figure('Position', [1 1 800 800])
% type = 3;
% hold on
% line([0 100],[0 100], 'color', 'k')
% for nf = 1:Nsgr
%     [~,i] = min(abs(long-coor_sugar(nf,2)));
%     [~,j] = min(abs(lat-coor_sugar(nf,1)));
%     scatter(para_sugar(type,nf), para_wrf(i,j,type), 150, 'fill')    
% end
% hold off
% grid on
% box  on
% 
% if(type == 1)
% xlim([0 3])
% ylim([0 3])
% xlabel('SuGAr amplitude (mm)')
% ylabel('WRF amplitude (mm)')
% elseif(type==2)
% set(gca,'xtick',0:3:24)
% set(gca,'ytick',0:3:24)
% xlim([0 24])
% ylim([0 24])
% xlabel('SuGAr phase (UTC)')
% ylabel('WRF phase (UTC)')
% elseif(type==3)
% xlim([46 56])
% ylim([46 56])
% set(gca,'xtick',46:2:58)
% set(gca,'ytick',46:2:58)
% xlabel('SuGAr mean (mm)')
% ylabel('WRF mean (mm)')
% end



