

load('../data/wrf_hgt.mat')
load('../data/wrf_domain.mat')
load('../data/gps_data_2008-2013_alph.mat')


coor_wx = [-0.875	100.352     2.7;
        1.166	97.705	6.1;
        1.556	98.889	13.1
        -3.864	102.339	15.2
        -5.242	105.179	88
        4.25	96.117	90];

N_wx    = size(coor_wx,1);    
   



gps_in   = importdata('../data/gps_ish_stations.txt');
N_gps    = size(gps_in.data,1);
coor_gps = zeros(N_gps,3);



for nf = 1:N_gps
    coor = cell2mat(gps_data(nf,2));
    coor_gps(nf,1) = coor(1);
    coor_gps(nf,2) = coor(2);
    
    [~,tmp] = min(abs( gps_in.data(nf)-coor_wx(:,3) ));
    
    coor_gps(nf,3) = tmp; %gps_in.data(nf);
end


figure('Position', [1 1 650 650])
hold on
contour(long,lat,hgt', 1:500:4001, 'linewidth', 1, 'color','k')

for nf = 1:N_wx
    scatter(coor_wx(nf,2),coor_wx(nf,1),1000,nf,  's', 'fill', 'markeredgecolor', 'k');
end

for nf = 1:N_gps
    scatter(coor_gps(nf,2),coor_gps(nf,1),250,coor_gps(nf,3),'fill', 'markeredgecolor', 'k');
end

hold off
box on
xlim([94 106]) %([94 110])
ylim([-6 6]) %([-8 8])
xlabel('Longitude')
ylabel('Latitude')
caxis([1 6])
colormap('parula')
