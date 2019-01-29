
load('wrf_domain.mat')

root  = '/n/regal/kuang_lab/gtorri/';
fn    = dir([root,'wrfout_d02*']);
N_wrf = length(fn);

Nz    = 44;

qrain = zeros(Nlat,Nlong,48);
ndays = 0;

for nf = 1:N_wrf
    
    ncp   = netcdf.open([root,fn(nf).name],'NOWRITE');
    
    ti = 1;
    if(nf == 1)
        ti = 481;
    end
    
    tf = 1440;
    Nt = 1440;
    if(nf == 13)
        tf = 277;
        Nt = 277;
    end    
    
    varid = netcdf.inqVarID(ncp,'QRAIN');
    tmp   = netcdf.getVar(ncp, varid, [0 0 0 0], [Nlat Nlong 1 Nt]);
    
    for t = ti:tf
        ts = rem(t-1,48)+1;
        qrain(:,:,ts)   = qrain(:,:,ts)  + tmp(:,:,t);
        disp(t)
    end
    ndays = ndays + (tf-ti+1)/48;
    disp(nf)
    
end

qrain = qrain/ndays;

save('wrf_qrain_avg.mat', 'qrain')





%%% Make movie

set(0, 'defaultaxesfontsize', 19)
set(0, 'defaultlinelinewidth', 2)


load('wrf_qrain_avg.mat')
load('wrf_domain.mat')



figure('Position', [1 1 650 650])
ts = 12;
hold on
contourf(long,lat,squeeze(qrain(:,:,ts)*1000)',200, 'linestyle', 'none')
contour(long,lat,lmask', 'k')
hold off
box on
xlim([95 106]) %([94 110])
ylim([-6 6]) %([-8 8])
caxis([0 0.25])
xlabel('Longitude')
ylabel('Latitude')

cb = colorbar;
title(cb, '(g kg^{-1})')

% mat = parula;
% jet_wrap = vertcat(flipud(mat(1:32,:)),flipud(mat(33:64,:)));
% colormap(jet_wrap);
%



% Cut temporal slice
loops = 48;
F(loops) = struct('cdata',[],'colormap',[]);
vidfile = VideoWriter('qrain_cycle.mp4','MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);

for ts = [(loops-14):loops 1:(loops-15)] %1:loops
    
    im = figure('Position', [1 1 650 650]);
    hold on
    contourf(long,lat,squeeze(qrain(:,:,ts)*1000)',200, 'linestyle', 'none')
    contour(long,lat,lmask', 'k')
    
    hold off
    
    xlim([95 106]) %([94 110])
    ylim([-6 6]) %([-8 8])
    caxis([0 0.25])
    
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    
    cb = colorbar;
    title(cb, '(g kg^{-1})')
    
    title(['Time LST = ', num2str(mod((ts+14)/2,24),'%02.1f'), ' hr'])
    
    drawnow
    F(ts) = getframe(gcf);
    
    writeVideo(vidfile,F(ts));
    
end
close(vidfile)
close all





