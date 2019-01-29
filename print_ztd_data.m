%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Accumulate data and compute averages  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear all
set(0, 'defaultaxesfontsize', 19)
set(0, 'defaultlinelinewidth', 2)


load('../data/gps_data_2008-2013_alph.mat')
load('../data/ztd_data_2009-2013_alph.mat')


load('../data/wrf_hgt.mat')
load('../data/wrf_domain.mat')
load('../data/distance_from_sumatra.mat')


Nf         = size(ztd_data,1);
coor_gps   = zeros(Nf,3);

IJ_gps     = zeros(Nf,2);
dist_gps   = zeros(Nf,1);

ztd_avg_dc = zeros(Nf,288);
num_avg_dc = zeros(Nf,288);


for nf = 1:Nf
    
    
    name = cell2mat(ztd_data(nf,1));
    for mf = 1:Nf
        name_tmp = cell2mat(gps_data(mf,1));
        if(name == name_tmp)
            
            coor = cell2mat(gps_data(nf,2));
            coor_gps(nf,1) = coor(1);
            coor_gps(nf,2) = coor(2);
            coor_gps(nf,3) = coor(3);
            
            [~,IJ_gps(nf,1)] = min( abs(coor(1) - lat) );
            [~,IJ_gps(nf,2)] = min( abs(coor(2) - long) );
            
            
            dist_gps(nf) = dist_mat(IJ_gps(nf,2),IJ_gps(nf,1));
            
            data = cell2mat(ztd_data(nf,2));
            
%             idx  = find(data(:,1) > 2009 & data(:,1) < 2013);
            idx  = find(data(:,1) == 2011 | data(:,1) == 2012);
            Nd   = length(idx);
            hr   = data(:,5);
            mn   = data(:,6);
            
            for m = 1:Nd
                nd = idx(m);
                
%                 if(data(nd,7) > 2350 && data(nd,7) < 2725)
                if(data(nd,7) > 2000)
                    
                    
                    tmp = (hr(nd))*12 + (mn(nd)/5) + 1;
                    
                    ztd_avg_dc(nf,tmp) = ztd_avg_dc(nf,tmp) + data(nd,7);
                    num_avg_dc(nf,tmp) = num_avg_dc(nf,tmp) + 1;
                end
            end
            
            disp(name)
            disp(nf)
            break
        end
    end
    
    
    
end


num_avg_dc(num_avg_dc < 11) = NaN;
ztd_avg_dc = ztd_avg_dc./num_avg_dc;


trsct = [37,35,27,26,1,34,32];
Nstat = length(trsct);


cmp = colormap(lines);

figure('Position', [1 1 700 650])
hold on
for m = 1:4
    vec = ztd_avg_dc(trsct(m),:);
    plot((1:288)/12, circshift( vec-mean(vec), 0) )
end
% for m = 5:Nstat
%     vec = ztd_avg_dc(trsct(m),:);
%     plot((1:288)/12, circshift( vec-mean(vec), 0), 'linestyle', '--')
% end


hold off

grid on 
box on

xlim([0 24])
set(gca, 'xtick', 0:3:24)

ylabel('ZTD anomaly (cm)')
xlabel('Time UTC (hr)')


% legend('SDKL','RNDG','PBLI','PBKR','ABGS','PTLO','PSMK', 'location', 'northwest')
legend('SDKL','RNDG','PBLI','PBKR', 'location', 'northwest')



% % Print all the diurnal cycles together
% Ncol = 22;
% cmp  = colormap(parula(Ncol));
% 
% figure('Position', [1 1 650 650])
% %for nf = 1:Nf
% for m = 1:Nstat
%     nf = trsct(m);
%     
%     hold on
%     
%     vec = ztd_avg_dc(nf,:);
%     dd  = floor(dist_gps(nf)/10) + 10;
%     
%     if(dd > 0 && dd < Ncol)
%         col = cmp(dd,:);
%         plot((1:288)/12, ( vec-mean(vec) ), 'color', col)
%     end
%     
% end
% hold off
% 
% box on
% grid on
% 
% ylabel('ZTD (cm)')
% xlabel('Time UTC (hr)')
% xlim([0 24])
% % ylim([-30 30])
% 
% 
% cb = colorbar;
% title(cb, '(km)')
% set(cb, 'xtick',0:0.1:1)
% set(cb, 'xticklabel',(0:0.1:1)*220-100)
% 
% cmp  = colormap(parula(Ncol));
% 
% 










