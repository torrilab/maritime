%%%%%%%%%%%%%%%%%%%%%%
%  A look using WRF  %
%%%%%%%%%%%%%%%%%%%%%%


set(0, 'defaultaxesfontsize', 19)
set(0, 'defaultlinelinewidth', 2)


rootz = '/n/regal/kuang_lab/gtorri/';

load([rootz,'wrf_hgt.mat']);
load([rootz,'wrf_domain.mat'])

min_lat  = min(lat);
max_lat  = max(lat);

min_long = min(long);
max_long = max(long);


dcycle_wrf = zeros(Nlong,Nlat,48);
ndays = 0;




U10    = zeros(Nlong,Nlat,48);
V10    = zeros(Nlong,Nlat,48);
QV_sfc = zeros(Nlong,Nlat,48);
TH_sfc = zeros(Nlong,Nlat,48);
P_sfc  = zeros(Nlong,Nlat,48);


for year = 2011:2012
    for nf = 1:13
        
        load([rootz,'/WRFOUT/',num2str(year),'/wrf_sfc_data_',num2str(nf),'.mat'])
        
        U10 = U10 + u10/13;
        V10 = V10 + v10/13;
        QV_sfc = QV_sfc + qv_sfc/13;
        TH_sfc = TH_sfc + th_sfc/13;
        P_sfc  = P_sfc  + p_sfc/13;
        disp(nf)
        
    end
end



U10_i  = zeros(Nlong+1,Nlat,48);
V10_i  = zeros(Nlong,Nlat+1,48);

for t = 1:48
    for j = 1:Nlat
        U10_i(:,j,t) = interp1(long,squeeze(U10(:,j,t)),long_i);
    end
    
    for i = 1:Nlong
        V10_i(i,:,t) = interp1(lat,squeeze(V10(i,:,t)),lat_i);
    end
end

conv = ((U10_i(1:Nlong,:,:)-U10_i(2:Nlong+1,:,:)) + (V10_i(:,1:Nlat,:)-V10_i(:,2:Nlat+1,:)))/3000;




for year = 2011:2012
    for nf = 1:13
        
        
        %     load(['../data/pwv_wrf_',num2str(nf),'_new.mat'])
        
        load([rootz,'/WRFOUT/',num2str(year),'/pwv_wrf_',num2str(nf),'_new.mat'])
        
        ti = 1;
        if(nf == 1)
            ti = 481;
        end
        
        tf = 1440;
        if(nf == 13)
            if(year == 2011)
                tf = 229;
            end
            if(year == 2012)
                tf = 277;
            end
        end
        
        for t = ti:tf
            ts = rem(t-1,48)+1;
            dcycle_wrf(:,:,ts) = dcycle_wrf(:,:,ts)+ pwv(:,:,t);
        end
        
        ndays = ndays + (tf-ti+1)/48;
        clear pwv;
        
        disp(nf)
        
    end
end

dcycle_wrf = dcycle_wrf/ndays;




% Select stations


% Transect 1
% MREK   2.969959559    98.538093958
% SDKL   2.796642489    98.192918679
% RNDG   2.665253991    97.857155506
% PBLI   2.308534775    97.405272011
% PBKR   2.065997053    97.082840126
% BABI   2.116813549    96.671201394

Nst = 6;

mask_1 = zeros(Nlat,Nlong);

trsc_1 = zeros(6,2);
trsc_1(1,:) = [2.969959559    98.538093958];
trsc_1(2,:) = [2.796642489    98.192918679];
trsc_1(3,:) = [2.665253991    97.857155506];
trsc_1(4,:) = [2.308534775    97.405272011];
trsc_1(5,:) = [2.065997053    97.082840126];
trsc_1(6,:) = [2.116813549    96.671201394];

for nf = 1:Nst
    vec = trsc_1(nf,:);
    
    [~,tmp_1]    = min(abs( vec(1)-lat  ));
    [~,tmp_2]    = min(abs( vec(2)-long ));
    trsc_1(nf,:) = [tmp_1 tmp_2];
end


for m = 1:Nst-1
    
    v_1 = trsc_1(m,:);
    v_2 = trsc_1(m+1,:);
    
    A = (v_2(1)-v_1(1))/(v_2(2)-v_1(2));
    
    for i = v_1(2):-1:v_2(2)
        j = floor(A*(i-v_1(2)) + v_1(1));
        mask_1(i,j) = 1;
    end
    
end


Nj = trsc_1(1,2) - trsc_1(Nst,2)+1;
pwv_avg  = zeros(Nj,48);
pwv_anm  = zeros(Nj,48);

u_avg    = zeros(Nj,48);
u_anm    = zeros(Nj,48);

v_avg    = zeros(Nj,48);
v_anm    = zeros(Nj,48);

conv_t   = zeros(Nj,48);

height   = zeros(1,Nj);

for j = trsc_1(Nst,2):trsc_1(1,2)
    
    idx   = find(mask_1(j,:) == 1);
    
    vec   = squeeze( dcycle_wrf(j,idx,:) );
    
    pwv_avg(j-trsc_1(Nst,2)+1,:) = vec';
    pwv_anm(j-trsc_1(Nst,2)+1,:) = (vec'-mean(vec));
    
    
    vec   = squeeze( U10(j,idx,:) );
    u_avg(j-trsc_1(Nst,2)+1,:) = vec';
    u_anm(j-trsc_1(Nst,2)+1,:) = (vec'-mean(vec));
    
    
    vec   = squeeze( V10(j,idx,:) );
    v_avg(j-trsc_1(Nst,2)+1,:) = vec';
    v_anm(j-trsc_1(Nst,2)+1,:) = (vec'-mean(vec));
    
    
    vec   = squeeze( conv(j,idx,:) );
    conv_t(j-trsc_1(Nst,2)+1,:) = vec';
    
    height(j-trsc_1(Nst,2)+1)    = hgt(j,idx);
    
end



cmp = colormap(lines(7));

figure('Position', [1 1 800 800])
hold on

% contour(long(trsc_1(Nst,2):trsc_1(1,2)),(0:47)/2,v_anm', -10:0.25:10)
contour(long(trsc_1(Nst,2):trsc_1(1,2)),(0:47)/2,pwv_anm', -4:0.5:4)
plot(long(trsc_1(Nst,2):trsc_1(1,2)),height/100, 'k')
for nf = 2:5
    scatter(long(trsc_1(nf,2)),    0, 200, 'fill', 'markeredgecolor','k', 'markerfacecolor', cmp(nf-1,:))
    scatter(long(trsc_1(nf,2)), 23.5, 200, 'fill', 'markeredgecolor','k', 'markerfacecolor', cmp(nf-1,:))
    line([long(trsc_1(nf,2)) long(trsc_1(nf,2))], [0 23.5], 'color', cmp(nf-1,:))
end

hold off
xlabel('Longitude')
ylabel('Time UTC (hr)')

colormap('redblue')
caxis([-4 4])
% caxis([-1 1])
box on
set(gca, 'ytick', 0:3:24)

cb = colorbar;
title(cb, '(mm)')




dcycle_anm_wrf = zeros(size(dcycle_wrf));
for i = 1:Nlong
    for j = 1:Nlat
        dcycle_anm_wrf(i,j,:) =  dcycle_wrf(i,j,:) - mean( dcycle_wrf(i,j,:),3);
    end
end


% Cut temporal slice
loops = 48;
F(loops) = struct('cdata',[],'colormap',[]);
vidfile = VideoWriter('movie_transect_1b.avi');
vidfile.FrameRate = 10;
open(vidfile);

for ts = 1:48
    im = figure('Position', [1 1 650 650]);
    
    hold on
    contour(long,lat, hgt', 1:100:1500, 'k')
%     contour(long,lat, squeeze(dcycle_anm_wrf(:,:,ts))', 100, 'linestyle', 'none')
    
    contour(long,lat, squeeze(dcycle_anm_wrf(:,:,ts))', -4:0.125:4)
    contour(long,lat, squeeze(dcycle_anm_wrf(:,:,ts))', -4:1:4, 'linewidth', 4)

    for m = 2:Nst-1
        scatter(long(trsc_1(m,2)),lat(trsc_1(m,1)), 200, 'markerfacecolor',cmp(m-1,:), 'markeredgecolor', 'k')
    end
    hold off
    box on
    
    ylim([1 3])
    xlim([96.5 98.5])
    %     caxis([35 55])
    caxis([-4 4])
    
    ylabel('Latitude (deg)')
    xlabel('Longitude (deg)')
    
    colormap(redblue)
    
    cb = colorbar;
    title(cb, '(mm)')
    
    title(['Time UTC: ', num2str((ts-1)/2), ' hr'])
    
    drawnow
    F(ts) = getframe(gcf);
    
    writeVideo(vidfile,F(ts));
    
end
close(vidfile)
close all








% Transect 2
% ABGS   0.220835382    99.387519481
% PTLO  -0.054589071    98.280033340
% PSMK  -0.089312892    97.860904240

mask_2 = zeros(Nlat,Nlong);


trsc_2 = zeros(3,2);
trsc_2(1,:) = [ 0.220835382    99.387519481];
trsc_2(2,:) = [-0.054589071    98.280033340];
trsc_2(3,:) = [-0.089312892    97.860904240];

for nf = 1:3
    vec = trsc_2(nf,:);
    
    [~,tmp_1]    = min(abs( vec(1)-lat  ));
    [~,tmp_2]    = min(abs( vec(2)-long ));
    trsc_2(nf,:) = [tmp_1 tmp_2];
end


for m = 1:2
    
    v_1 = trsc_2(m,:);
    v_2 = trsc_2(m+1,:);
    
    A = (v_2(1)-v_1(1))/(v_2(2)-v_1(2));
    
    for i = v_1(2):-1:v_2(2)
        j = floor(A*(i-v_1(2)) + v_1(1));
        mask_2(i,j) = 1;
    end
    
end






% figure('Position', [1 1 800 800])
% hold on
% % contourf(long,lat,mask_1'+mask_2',200, 'linestyle', 'none')
% pcolor(long,lat,mask_1'+mask_2')
%
% contour(long,lat,hgt', 1000:1000:4000, 'linewidth', 1, 'color','r')
% contour(long,lat,hgt', 1:60:300, 'linewidth', 1, 'color','k')
% hold off
% box on
% xlim([95 106]) %([94 110])
% ylim([-6 6]) %([-8 8])
% xlabel('Longitude')
% ylabel('Latitude')
% caxis([0 1])
% cb = colorbar;
%






