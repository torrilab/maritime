clear all

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultaxesfontsize', 19)

load('../data/wrf_domain.mat')


% Find the connected components of the figures and select the biggest
CC            = bwconncomp(lmask);
numPixels     = cellfun(@numel,CC.PixelIdxList);
[~,idx]       = max(numPixels);
Sumatra       = CC.PixelIdxList{idx};
[Sum_i,Sum_j] = ind2sub([Nlong,Nlat],Sumatra);


% Draw the island of sumatra
s_mask = zeros(Nlong,Nlat);
t_mask = zeros(Nlong,Nlat);
u_mask = zeros(Nlong,Nlat);
for n = 1:length(Sum_i)
    i = Sum_i(n);
    j = Sum_j(n);
    
    s_mask(i,j) = 1;
end


% Draw its contour
B       = bwboundaries(s_mask);
mat_bnd = B{1};
for n = 1:size(mat_bnd,1)
    i = mat_bnd(n,1);
    j = mat_bnd(n,2);
    
    t_mask(i,j) = 1;
end


% Select the bottom-right and top-left points
val_x = min(mat_bnd(:,2));
idi   = find(mat_bnd(:,2) == val_x);
val_y = max(mat_bnd(idi,1));
BR    = [val_y,val_x];

val_y = min(mat_bnd(:,1));
idj   = find(mat_bnd(:,1) == val_y);
val_x = max(mat_bnd(idj,2));
TL    = [val_y,val_x];


% Cut out the two extermes
for a1 = -1:2
    for a2 = -1:1
        i = TL(1);
        j = TL(2);
        t_mask(i+a1,j+a2) = 0;
        
        i = BR(1);
        j = BR(2);
        t_mask(i+a1,j+a2) = 0;
    end
end


% Find the smallest connected component of the boundary
DD            = bwconncomp(t_mask);
numPixels     = cellfun(@numel,DD.PixelIdxList);
[~,idx]       = min(numPixels);
coast         = DD.PixelIdxList{idx};
[cst_i,cst_j] = ind2sub([Nlong,Nlat],coast);
for n = 1:length(coast)
    i = cst_i(n);
    j = cst_j(n);
    u_mask(i,j) = 1;
end
dist_mat = bwdist(u_mask)*3;


% Make distances negative to the norht-east of the island
for j = 1:Nlat
    for i = 1:Nlong
        
        X = [i,j];
        Y = [cst_i,cst_j];
        ind = dsearchn(Y, X);
        tmp = Y(ind,:)-X;
        if(tmp(1) < 0 || tmp(2) < 0)
            dist_mat(i,j) = -dist_mat(i,j);
        end
        
    end
end


% Read WRF input
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


% Order things based on distance from the coast
max_dist = 500;
min_dist = -500;
d_dist   = 10;
bins     = min_dist:d_dist:max_dist;
Nbins    = length(bins);
pwv_dist = zeros(Nbins,48);
numt     = zeros(1,Nbins);

tmp = dist_mat;
tmp(tmp < min_dist) = NaN;
tmp(tmp > max_dist) = NaN;
tmp1 = floor((tmp - min_dist)/d_dist)+1;

for j = 1:Nlat
    for i = 1:Nlong
        pt = tmp1(i,j);
        if(~isnan(pt))
            pwv_dist(pt,:) = pwv_dist(pt,:) + squeeze( dcycle_wrf(i,j,:) )';
            numt(pt) = numt(pt) + 1;
        end
    end
end

pwv_dist_avg = pwv_dist./numt';



figure('Position', [1 1 800 800])
hold on
contourf(-bins(1:end),(0:47)/2,circshift(pwv_dist_avg,[0 14])', 38:0.2:55, 'linestyle', 'none')
contour(-bins(1:end),(0:47)/2,circshift(pwv_dist_avg,[0 14])', 47:0.2:62, 'k')
hold off
box on

xlim([-100 20])
set(gca, 'xtick', -100:20:20)
set(gca, 'xticklabel', 100:-20:-20)

% xlim([-280 20])
% set(gca, 'xtick', -280:50:20)
% set(gca, 'xticklabel', 280:-50:-20)

set(gca, 'ytick', 0:3:24)

caxis([45 51])
cb = colorbar;
title(cb, '(mm)')

xlabel('Distance from coast (km)')
ylabel('LST (hr)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Close-up on a particular area of Sumatra %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% lat  range: -2 -> -0.8
% long range: 98 -> 100

[~,j_i] = min(abs(lat - (-2)));
[~,j_f] = min(abs(lat - (-0.8)));

[~,i_i] = min(abs(long - (98)));
[~,i_f] = min(abs(long - (100)));

pwv_zoom = zeros(2,48);
num_zoom = zeros(2,48);

for j = j_i:j_f
    for i = i_i:i_f

        c = lmask(i,j)+1;
        pwv_zoom(c,:) = pwv_zoom(c,:) + squeeze( dcycle_wrf(i,j,:) )';
        num_zoom(c,:) = num_zoom(c,:) + 1;

    end
end
pwv_zoom = pwv_zoom./num_zoom;



% Compare diurnal cycles, land vs. ocean
figure('Position', [1 1 800 800])
plot(0:0.5:23.5,circshift(pwv_zoom, [0 14])')
grid on 
box on
xlim([0 24])
set(gca, 'xtick', 0:3:24)
xlabel('LST (hr)')
ylabel('PWV (mm)')
legend('Ocean','Land')






% Cut longitudinal slice 
i_tmp = floor((i_f+i_i)/2);
mat_tmp = squeeze( dcycle_wrf(i_tmp,:,:) );

% vec_tmp = lmask(i_tmp,j_i:j_f);
% [~,pos_i] = max( vec_tmp(2:end)-vec_tmp(1:end-1) );
% [~,pos_f] = min( vec_tmp(2:end)-vec_tmp(1:end-1) );
% pos_i = pos_i + 1;

figure('Position', [1 1 800 800])
hold on
contourf(0:0.5:23.5,lat,circshift(mat_tmp,[0,14]), 46:0.2:54)
% line([-1000 1000],[lat(pos_i+j_i-1) lat(pos_i+j_i-1)],  'color','r')
% line([-1000 1000],[lat(pos_f+j_i-1) lat(pos_f+j_i-1)],  'color','r')
plot(hgt(i_tmp,j_i:j_f)/100,lat(j_i:j_f), 'color', 'k')
hold off
ylim([-2 -0.8])
xlim([0 23.5])

set(gca, 'xtick', 0:3:24)
xlabel('LST (hr)')
ylabel('Latitude (deg)')

cb = colorbar;
title(cb, '(mm)')


% Cut latitudinal slice 
j_tmp = floor((j_f+j_i)/2);
mat_tmp = squeeze( dcycle_wrf(:,j_tmp,:) );

% vec_tmp = lmask(i_i:i_f,j_tmp);
% [~,pos_i] = max( vec_tmp(2:end)-vec_tmp(1:end-1) );
% [~,pos_f] = min( vec_tmp(2:end)-vec_tmp(1:end-1) );
% pos_i = pos_i + 1;

figure('Position', [1 1 800 800])
hold on
contourf(0:0.5:23.5,long,circshift(mat_tmp,[0,14]), 46:0.2:54)
plot(hgt(:,j_tmp)/100,long, 'color', 'black')
hold off
ylim([98 101])
xlim([0 23.5])
caxis([46 54])
set(gca, 'xtick', 0:3:24)
xlabel('LST (hr)')
ylabel('Longitude (deg)')

cb = colorbar;
title(cb, '(mm)')

view(90,-90)




% Show mask with slices cut
figure('Position', [1 1 800 800])
hold on
contourf(long,lat,lmask')
line([long(i_tmp) long(i_tmp)], [-1000 1000], 'color','r')
line([-1000 1000], [lat(j_tmp) lat(j_tmp)], 'color','g')
hold off
xlim([98 100])
ylim([-2 -0.8])
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')



% Cut temporal slice
loops = 48;
F(loops) = struct('cdata',[],'colormap',[]);
vidfile = VideoWriter('siberut_large.mp4','MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);

for ts = [(loops-14):loops 1:(loops-15)]
    
    mat_tmp = squeeze( dcycle_wrf(:,:,ts) );
%     im = figure('Position', [1 1 800 800]);
    im = figure('Position', [1 1 1200 600]);
    hold on
    contourf(long,lat,mat_tmp',  46:0.25:54, 'linestyle', 'none')
    contour(long,lat, hgt', 5:50:1000, 'k', 'linewidth', 1)
% 
%     contourf(long,lat,mat_tmp',  0:0.25:54, 'linestyle', 'none')
%     contour(long,lat, hgt', 5:250:10000, 'k', 'linewidth', 1)

    hold off
    
    xlim([97 101])
    ylim([-2 0])
    caxis([46 53])
    
%     xlim([95 106])
%     ylim([-6 6])
%     caxis([30 54])

    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    
    cb = colorbar;
    title(cb, '(mm)')
    
    title(['Time LST = ', num2str(mod((ts+14)/2,24),'%02.1f'), ' hr'])
    
    drawnow
    F(ts) = getframe(gcf);
    
    writeVideo(vidfile,F(ts));
    
end
close(vidfile)
close all





