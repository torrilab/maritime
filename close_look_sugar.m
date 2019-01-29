%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Accumulate data and compute averages  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear all
load('gps_data_2010-2012.mat')


Nf   = size(gps_data,1);
nf   = Nf-6; %12

name = cell2mat(gps_data(nf,1));
coor = cell2mat(gps_data(nf,2));
data = cell2mat(gps_data(nf,3));



hrs  = data(:,5);
pwv  = data(:,8);


Ndata = size(pwv,1);

Ndays = -999;
for i = 1:8
    tmp   = find(hrs == 3*(i-1));
    disp(length(tmp))
    Ndays = max(Ndays,length(tmp));
end


% Divide the data into days
pwv_days = zeros(Ndays,8);

nd  = 1;
day = 1;
kr  = floor(hrs(nd)/3)+1;
pwv_days(day,kr) = pwv(nd);
for nd = 2:Ndata
    if(hrs(nd) < hrs(nd-1))
        day = day + 1;
    end
    kr  = floor(hrs(nd)/3)+1;
    pwv_days(day,kr) = pwv(nd);
end




% Find maximum
Ndays = size(pwv_days,1);
histo_max = zeros(1,8);
for nd = 1:Ndays
    
    vec = pwv_days(nd,:);
    id0 = find(vec == 0);
    
    if(length(id0) < 2)
        [val,pos] = max(vec);
        histo_max(pos) = histo_max(pos)+1;
    end
    
end



figure('Position', [1 1 800 800])
hold on
b1 = bar((0:3:21), [ (lais_histo_max)/sum(lais_histo_max) ; (ngng_histo_max)/sum(ngng_histo_max) ]' );
% b2 = bar((0:3:21), ngng_histo_max);
hold off
grid on
box on
xlim([-1 22])

b1(1).FaceColor =  [0    0.4470    0.7410];
b1(2).FaceColor = [0.8500    0.3250    0.0980];
legend('LAIS', 'NGNG')

set(gca,'xTick', 0:3:21)

xlabel('Time UTC (hrs)')
ylabel('Fraction of days')

