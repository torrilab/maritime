clear all
set(0, 'DefaultAxesFontSize', 19)


% Read GPS stations
nm      = 'PSMK';
yr      = 2013;

coord   = [-0.089312892    97.860904240      13.3302868];
dowrite = 0;
domap   = 0;
coeff   = -1;
% name    = ['data_',num2str(yr),'/data_wx_',nm,'_',num2str(yr),'.txt'];
name    = ['data_wx_',nm,'_',num2str(yr),'_new.txt'];

% Read weather stations
data_wx = xlsread('../data/ish/station-Sumatra.xlsx');
N_raw   = size(data_wx,1);
w_raw   = size(data_wx,2);
wx_stat = zeros(0,w_raw);

for n = 1:N_raw
    
    y_ini = floor(data_wx(n,9)/10000);
    y_fin = floor(data_wx(n,10)/10000);
    
    % Select only those stations that include time of interest
    if(y_ini <= 2010 && y_fin >= 2012)
            wx_stat = cat(1,wx_stat,data_wx(n,:));
    end
    
end

N_wx   = size(wx_stat,1);
wx_coo = wx_stat(:,6:8);


%%% WX %%%

% Rank weather stations according to distance
wx_dist = sqrt((wx_coo(:,1)-coord(1)).^2 + (wx_coo(:,2)-coord(2)).^2);


% Find closest weather station
[val,pos] = min(abs(wx_dist));
pos=pos+coeff
data_str  = [num2str(wx_stat(pos,1)),'-', num2str(wx_stat(pos,2)),'-', num2str(yr)];


% Collect data from that station
wx_raw  = importdata(['../data/ish/Sumatra/',data_str]);
Nd_wx   = length(wx_raw);
data_wx = zeros(Nd_wx,5);

for n = 1:Nd_wx
    
    data_tmp1   = cell2mat(wx_raw(n));
    
    vec_year    = str2double(data_tmp1(16:19));
    vec_month   = str2double(data_tmp1(20:21));
    vec_day     = str2double(data_tmp1(22:23));
    vec_hour    = str2double(data_tmp1(24:25));
    vec_minute  = str2double(data_tmp1(26:27));
    vec_second  = 00;
    
    data_wx(n,1) = datenum([ vec_year, vec_month, vec_day, ...
        vec_hour, vec_minute, vec_second]);        % Time
    
    data_wx(n,2) = str2double(data_tmp1(88:92));   % Temperature (C)
    data_wx(n,3) = str2double(data_tmp1(93));      % Quality index
    
    data_wx(n,4) = str2double(data_tmp1(100:104)); % Pressure (hPa)
    data_wx(n,5) = str2double(data_tmp1(105));     % Quality index
    
end

data_wx(data_wx == 99999) = NaN;
data_wx(data_wx == 9999) = NaN;
data_wx(:,1) = (data_wx(:,1) - data_wx(1,1))*24;
data_wx(:,2) = data_wx(:,2)/10;
data_wx(:,4) = data_wx(:,4)/10;

disp(wx_coo(pos,:))

if(dowrite)
dlmwrite(name, data_wx)
end
  
if(domap)
% Plot all the available stations that cover 2010 and 2012
figure('Position', [1 1 800 800])
hold on
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
scatter(wx_coo(:,2),wx_coo(:,1),    150, 'fill')
scatter(coord(:,2),coord(:,1),       150, 'fill')
scatter(wx_coo(pos,2),wx_coo(pos,1), 150)
grid on
box on
xlim([94 110])
ylim([-8 8])
xlabel('Longitude')
ylabel('Latitude')


num   = zeros(1,8);
T_avg = zeros(1,8);
p_avg = zeros(1,8);
for nd = 1:Nd_wx

    ix = floor(mod(data_wx(nd,1),24)/3 + 1);
    
    if(~isnan(data_wx(nd,2)) && ~isnan(data_wx(nd,4)))
        num(ix) = num(ix) + 1;
        
        T_avg(ix) = T_avg(ix) + data_wx(nd,2);
        p_avg(ix) = p_avg(ix) + data_wx(nd,4);
    end
    
end
T_avg = T_avg./num;
p_avg = p_avg./num;




figure('Position', [1 1 800 600])

yyaxis left
plot((0:8)*3+1,cat(2,circshift(T_avg,2),T_avg(end-1)))
ylabel('Temperature (C)')

yyaxis right
plot((0:8)*3+1,cat(2,circshift(p_avg,2),p_avg(end-1)))
ylabel('Pressure (mb)')
xlim([1 25])

grid on
set(gca, 'xtick', 1:3:25)
set(gca, 'xticklabel', [1:3:22 1])
xlabel('LST (hr)')



end





