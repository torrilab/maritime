
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A look using the GPS stations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../data/gps_data_2008-2013_alph.mat')

% MREK - SDKL - RNDG - PBLI - PBKR - BABI
%  NaN -  37  -  35  -  27  -  26  -  NaN


% ABGS - LBHU - PBAI - PTLO - PSMK
%   1  -  NaN -  NaN -  34  -  32



trsct = [37,35,27,26,1,34,32];
Nstat = length(trsct);

gps_trsct = cell(Nstat,3);
for n = 1:Nstat
    gps_trsct(n,:) = gps_data(trsct(n),:);
end


dcycle_trsct = zeros(Nstat,8);
numb_trsct   = zeros(Nstat,8);

for nf = 1:Nstat
    
    data = cell2mat(gps_trsct(nf,3));
    
    yr   = data(:,1);
    hr   = data(:,5);
    
    idx  = find(data(:,1) == 2011 | data(:,1) == 2012);
    Nd   = length(idx);
%     for nd = 1:Nd
    for m = 1:Nd
        nd = idx(m);
        tmp = floor(hr(nd)/3)+1;
        dcycle_trsct(nf,tmp) = dcycle_trsct(nf,tmp) + data(nd,8);
        numb_trsct(nf,tmp) = numb_trsct(nf,tmp) + 1;
    end
    
end
numb_trsct(numb_trsct < 5) = NaN;
tmp_1 = nansum(dcycle_trsct,2);
tmp_2 = nansum(numb_trsct,2);
vec   = (tmp_1./tmp_2);

dcycle_trsct = dcycle_trsct./numb_trsct;


cmp = colormap(lines);

figure('Position', [1 1 700 650])
mat = cat(2, dcycle_trsct, dcycle_trsct(:,1));
hold on
for nf = 1:4
plot((0:8)*3, mat(nf,:)-mean(mat(nf,:)), 'linewidth', 2)
end

for nf = 5:Nstat
% plot((0:8)*3, mat(nf,:)-mean(mat(nf,:)), 'linestyle', '--')
end

hold off
grid on
box on 

xlim([0 24])
set(gca, 'xtick', 0:3:24)

ylabel('PWV anomaly (mm)')
xlabel('Time UTC (hr)')
% legend('SDKL','RNDG','PBLI','PBKR','ABGS','PTLO','PSMK', 'location', 'northwest')
legend('SDKL','RNDG','PBLI','PBKR', 'location', 'northwest')









