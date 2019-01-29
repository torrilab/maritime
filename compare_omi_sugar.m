
clear all


load('../data/gps_data_2008-2013_alph.mat')

sugar_name = cell2mat(gps_data(:,1));

N_sgr = size(sugar_name,1);

dr = dir('../data/OMI/2011/');

N_omi = length(dr);

% datum = zeros(0,2);

d_sugar = zeros(1,0);
d_omi   = zeros(1,0);

for nomi = 3:N_omi
    
    nfile  = dr(nomi).name;
    
    omi_in = importdata(['../data/OMI/2011/',nfile]);
    mat    = omi_in.data;
    
    Ndays  = size(mat,1);
    
    whn    = mat(:,1);
    wht    = mat(:,2);
    
    idx    = find(wht > 0);
    
    if(~isempty(idx))
        
        year   = floor(whn(idx)/10000);
        month  = floor((whn(idx) - year*10000)/100);
        day    = whn(idx) - year*10000 - month*100;
        
        dtmp1  = wht(idx);
        dtmp2  = size(dtmp1);
        
%         d_omi  = cat(2,dtmp1',d_omi);
        
        for nsgr = 1:N_sgr
            if(nfile(1:4) == lower(sugar_name(nsgr,:)))
                disp(nsgr)
                break
            end
        end
        
        sugar_data = cell2mat(gps_data(nsgr,3));
        for m = 1:length(idx)
            n = idx(m);
            
            vec = find( sugar_data(:,1) == year(m) & sugar_data(:,3) == month(m) & sugar_data(:,4) == day(m));
            if(~isempty(vec))
                dtmp2(m) = sugar_data(vec(2),8);
                
                d_omi   = cat(2,dtmp1(m),d_omi);
                d_sugar = cat(2,dtmp2(m),d_sugar);
                
%                 disp([nomi dtmp1(m) dtmp2(m)])
            end
        end
        
        clear dtmp1 dtmp2

    end
    
end





figure('Position', [1 1 650 650])
hold on
line([0 100], [0 100], 'color', [0.8500    0.3250    0.0980], 'linewidth', 2)
scatter(d_sugar,d_omi, 250,'fill', 'markeredgecolor', 'k')
hold off
grid on
box on

xlim([20 80])
ylim([20 80])

xlabel('PWV from SuGAr (mm)')
ylabel('PWV from OMI (mm)')





