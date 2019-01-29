
load('wrf_domain.mat')

Nlat  = 402;
Nlong = 402;
Nts   = 48;
N_wrf = 6;

u10_tot    = zeros(Nlat,Nlong,Nts);
v10_tot    = zeros(Nlat,Nlong,Nts);
p_sfc_tot  = zeros(Nlat,Nlong,Nts);
qv_sfc_tot = zeros(Nlat,Nlong,Nts);
th_sfc_tot = zeros(Nlat,Nlong,Nts);
sst_tot    = zeros(Nlat,Nlong,Nts);



for nf  = [1:4 6:13]

    load(['wrf_sfc_data_',num2str(nf),'.mat'])
    
    u10_tot    = u10_tot + u10/N_wrf;
    v10_tot    = v10_tot + v10/N_wrf;
    p_sfc_tot  = p_sfc_tot + p_sfc/N_wrf;
    qv_sfc_tot = qv_sfc_tot + qv_sfc/N_wrf;
    th_sfc_tot = th_sfc_tot + th_sfc/N_wrf;
    sst_tot    = sst_tot + sst/N_wrf;
    
end




% Cut latitudinal slice 
j_tmp = floor((j_f+j_i)/2);
mat_tmp = squeeze( u10_tot(:,j_tmp,:) );

% vec_tmp = lmask(i_i:i_f,j_tmp);
% [~,pos_i] = max( vec_tmp(2:end)-vec_tmp(1:end-1) );
% [~,pos_f] = min( vec_tmp(2:end)-vec_tmp(1:end-1) );
% pos_i = pos_i + 1;

figure('Position', [1 1 800 800])
hold on
contourf(0:0.5:23.5,long,mat_tmp,-7:0.25:7)
plot(hgt(:,j_tmp)/100,long, 'color', 'black')
hold off
ylim([98 101])
xlim([0 23.5])
caxis([-6 6])
set(gca, 'xtick', 0:3:24)
xlabel('Time UTC (hr)')
ylabel('Longitude (deg)')
colormap('redblue')
cb = colorbar;
title(cb, '(mm)')

view(90,-90)






% Cut temporal slice
loops = 48;
F(loops) = struct('cdata',[],'colormap',[]);
vidfile = VideoWriter('siberut_large_u.mp4','MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);

for ts = 1:loops
    
    mat_tmp = squeeze( u10(:,:,ts) );
%     im = figure('Position', [1 1 800 800]);
    im = figure('Position', [1 1 1200 600]);
    hold on
    contourf(long,lat,mat_tmp',  -7:0.25:7, 'linestyle', 'none')
    contour(long,lat, hgt', 5:50:1000, 'k', 'linewidth', 1)

    
    hold off
    
    xlim([97 101])
    ylim([-2 0])
    caxis([-6 6])
    
%     xlim([95 106])
%     ylim([-6 6])
%     caxis([-7 7])

    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    colormap('redblue')
    cb = colorbar;
    title(cb, '(m s^{-1})')
    
    title(['Time UTC = ', num2str(ts/2,'%02.1f'), ' hr'])
    
    drawnow
    F(ts) = getframe(gcf);
    
    writeVideo(vidfile,F(ts));
    
end
close(vidfile)
close all


