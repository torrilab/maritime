clear all

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultaxesfontsize', 19)

branch = '/n/regal/kuang_lab/gtorri/';

load([branch,'wrf_hgt.mat']);
load([branch,'wrf_domain.mat']);


dcycle_wrf = zeros(Nlat,Nlong,48);
ndays = 0;

for year = 2011 %:2012 
for nf = 1:13
    
    load([branch,'/WRFOUT/',num2str(year),'/pwv_wrf_',num2str(nf),'_new.mat'])
    
    ti = 1;
    if(nf == 1)
        ti = 481;
    end
    
    tf = 1440;
    if(nf == 13)
    if(year == 2011)
        tf = 229;
    end
    if(year == 2012)
        tf = 277;
    end
    end
    
    for t = ti:tf
        ts = rem(t-1,48)+1;
        dcycle_wrf(:,:,ts) = dcycle_wrf(:,:,ts)+ pwv(:,:,t);
    end
    
    ndays = ndays + (tf-ti+1)/48;
    clear pwv;
    
    disp(nf)
    
end
end

dcycle_wrf = dcycle_wrf/ndays;


para_wrf = zeros(Nlat,Nlong,3);
x_wrf    = 2*pi*(0:0.5:23.5)/24;
for j = 1:Nlong
    for i = 1:Nlat
        
        % Read the data in a format to feed the regression
        y   = squeeze( dcycle_wrf(i,j,:) )';
        ym  = mean(y);
        yr  = (max(y)-min(y));
        fit = @(b,x_wrf)  b(1).*(cos(x_wrf - b(2))) + b(3);
        fcn = @(b) sum((fit(b,x_wrf) - y).^2);
        para_wrf(i,j,:) = fminsearch(fcn, [yr;  0;  ym]);
        
    end
    
    disp(j)
end
para_wrf(:,:,2) = para_wrf(:,:,2)/(2*pi)*24 + 7*0;


% Impose that amplitude be positive
vec1 = para_wrf(:,:,1);
vec2 = para_wrf(:,:,2);

vec2(vec1 < 0) = vec2(vec1 < 0)+12;
vec1 = abs(vec1);
para_wrf(:,:,1) = vec1;

vec2(vec2 < 0)  = vec2(vec2 < 0)  + 24;
vec2(vec2 > 24) = vec2(vec2 > 24) - 24;
para_wrf(:,:,2) = vec2;





save('dcycle_wrf_2011_UTC.mat', 'dcycle_wrf', 'para_wrf')







