clear all

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultaxesfontsize', 19)

load('../data/wrf_hgt.mat')
load('../data/wrf_domain.mat')


dcycle_wrf = zeros(Nlat,Nlong,48);
ndays = 0;

for nf = 1:13
    
    load(['../data/pwv_wrf_',num2str(nf),'_new.mat'])
    
    ti = 1;
    if(nf == 1)
        ti = 481;
    end
    
    tf = 1440;
    if(nf == 13)
        tf = 277;
    end
    
    for t = ti:tf
        ts = rem(t-1,48)+1;
        dcycle_wrf(:,:,ts) = dcycle_wrf(:,:,ts)+ pwv(:,:,t);
    end
    
    ndays = ndays + (tf-ti+1)/48;
    clear pwv;
    
    disp(nf)
    
end
dcycle_wrf = dcycle_wrf/ndays;


para_wrf = zeros(Nlat,Nlong,3);
x_wrf    = 2*pi*(0:0.5:23.5)/24;
for j = 1:Nlong
    for i = 1:Nlat
        
        % Read the data in a format to feed the regression
        y   = squeeze( dcycle_wrf(i,j,:) )';
        ym  = mean(y);
        yr  = (max(y)-min(y));
        fit = @(b,x_wrf)  b(1).*(cos(x_wrf - b(2))) + b(3);
        fcn = @(b) sum((fit(b,x_wrf) - y).^2);
        para_wrf(i,j,:) = fminsearch(fcn, [yr;  0;  ym]);
        
    end
    
    disp(j)
end
para_wrf(:,:,2) = para_wrf(:,:,2)/(2*pi)*24 + 7;


% Impose that amplitude be positive
vec1 = para_wrf(:,:,1);
vec2 = para_wrf(:,:,2);

vec2(vec1 < 0) = vec2(vec1 < 0)+12;
vec1 = abs(vec1);
para_wrf(:,:,1) = vec1;

vec2(vec2 < 0)  = vec2(vec2 < 0)  + 24;
vec2(vec2 > 24) = vec2(vec2 > 24) - 24;
para_wrf(:,:,2) = vec2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute climatology for SuGAr %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('../data/gps_data_2010-2012.mat')

Nsgr = size(gps_data,1);

coor_sugar   = zeros(Nsgr,2);
dcycle_sugar = zeros(Nsgr,8);
num_sugar    = zeros(Nsgr,8);

for nf = 1:Nsgr
    
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
para_sugar(2,:) = para_sugar(2,:)/(2*pi)*24 + 7;


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
    
%     dcycle_point_wrf(:,m)   = dcycle_wrf(i,j,:);
    
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
        
    load(['../data/pwv_wrf_',num2str(ng),'_new.mat'])
    
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Print output figures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Print average of PWV
type = 2;

figure('Position', [1 1 800 800])
hold on
contourf(long,lat,squeeze(para_wrf(:,:,type))',200, 'linestyle', 'none')
for nf = 1:Nsgr
    scatter(coor_sugar(nf,2),coor_sugar(nf,1),200,para_sugar(type,nf), 'fill', 'markeredgecolor', 'k')
end
contour(long,lat,hgt', 1000:1000:4000, 'linewidth', 1, 'color','r')
contour(long,lat,hgt', 1:60:300, 'linewidth', 1, 'color','k')
hold off
box on
xlim([95 106]) %([94 110])
ylim([-6 6]) %([-8 8])
xlabel('Longitude')
ylabel('Latitude')

cb = colorbar;


% title(cb, '(mm)')
% caxis([0 4])

title(cb, '(LST)')
caxis([0 24])
set(cb, 'xtick', 0:3:24)

mat = parula;
jet_wrap = vertcat((mat(1:64,:)),flipud(mat(1:64,:)));
colormap(jet_wrap);





figure('Position', [1 1 800 800])
type = 2;
hold on
line([0 100],[0 100], 'color', 'k')
for nf = 1:Nsgr
    
    [~,i] = min(abs(long-coor_sugar(nf,2)));
    [~,j] = min(abs(lat-coor_sugar(nf,1)));
    
    scatter(para_sugar(type,nf), para_wrf(i,j,type), 150, 'fill')
    
end
hold off
grid on
box  on

if(type == 1)
xlim([0 3])
ylim([0 3])
xlabel('SuGAr amplitude (mm)')
ylabel('WRF amplitude (mm)')
elseif(type==2)
set(gca,'xtick',0:3:24)
set(gca,'ytick',0:3:24)
xlim([0 24])
ylim([0 24])
xlabel('SuGAr phase (LST)')
ylabel('WRF phase (LST)')
elseif(type==3)
xlim([46 56])
ylim([46 56])
set(gca,'xtick',46:2:58)
set(gca,'ytick',46:2:58)
xlabel('SuGAr mean (mm)')
ylabel('WRF mean (mm)')
end



