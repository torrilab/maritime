clear all

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultaxesfontsize', 19)


root = '/n/regal/kuang_lab/gtorri/WRFOUT/2012/';

load('../data/wrf_domain.mat')


% Read WRF input
dcycle_wrf = zeros(Nlat,Nlong,48);
ndays = 0;

for nf = 1:13
    
    load([root,'/pwv_wrf_',num2str(nf),'_new.mat'])
    
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




