

clear all
load('../data/dcycle_wrf_comp.mat')
load('../data/dcycle_nvap_comp.mat')


wrf_t   = 0:0.5:24;
sugar_t = 0:3:24;
nvap_t  = 0:6:24;


wrf_pwv_tmp   = cat(1, dcycle_wrf_comp,   dcycle_wrf_comp(1,:));
sugar_pwv_tmp = cat(1, dcycle_sugar_comp, dcycle_sugar_comp(1,:));
nvap_pwv_tmp  = cat(1, dcycle_nvap_comp,  dcycle_nvap_comp(1,:));


wrf_pwv   = wrf_pwv_tmp;
sugar_pwv = interp1(sugar_t, sugar_pwv_tmp, wrf_t);
nvap_pwv  = interp1(nvap_t,  nvap_pwv_tmp,  wrf_t);

color_m   = [[0         0.4470    0.7410]', ...
             [0.8500    0.3250    0.0980]', ...
             [0.9290    0.6940    0.1250]', ...
             [0.4940    0.1840    0.5560]']';


         
         
figure('Position', [1 1 1200 600])

sb1 = subplot(1,2,1);
hold on
for n = 1:4
plot((0:0.5:24),sugar_pwv(:,n), 'color', color_m(n,:))
end

for n = 1:4
plot((0:0.5:24),nvap_pwv(:,n),  'color', color_m(n,:), 'linestyle', '--')
end
hold off
grid on
box on
xlim([0 24])
ylim([28 60])
set(gca, 'xtick', 0:3:24)
xlabel('UTC (hrs)')
ylabel('PWV (mm)')
sb1.Position = [0.07 0.11 0.395 0.8150];

sb2 = subplot(1,2,2);
hold on
for n = 1:4
plot((0:0.5:24),sugar_pwv(:,n), 'color', color_m(n,:))
end

for n = 1:4
plot((0:0.5:24),wrf_pwv(:,n),   'color', color_m(n,:), 'linestyle', ':')
end
hold off
grid on
box on
xlim([0 24])
ylim([28 60])
set(gca, 'xtick', 0:3:24)
xlabel('UTC (hrs)')
ylabel('PWV (mm)')
legend('SDKL','RNDG','PBLI','PBKR', 'location', 'southwest')
sb2.Position = [0.55 0.11 0.395 0.8150];

