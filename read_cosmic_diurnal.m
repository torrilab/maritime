clear all
set(0, 'DefaultAxesFontSize', 19)
set(0, 'DefaultLineLineWidth', 2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define the GPS stations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gps_chosen = [
        2.520820653    96.150823673      14.2806945 %BNON
        2.409254108    96.326161460     -18.3090609 %BSIM
        2.923604066    95.804053163      12.9097333 %LEWK
    
        1.078615694    97.811361410      -3.1571225 %BITI
        0.569193443    97.710637087      73.1791457 %BTHL
    
       -3.076698795   100.284551202       8.2304011 %BSAT
       -3.529233353   102.033934771      23.0053792 %LAIS
       
        0.220835382    99.387519481     241.7700422 %ABGS
       -1.281540959    98.643940064      23.0071832 %BTET
        2.796642489    98.192918679     894.4091617 %SDKL
        
        2.665253991    97.857155506       3.3826798 %RNDG
       -2.966606371   100.399604988      23.5023307 %PBKR
        2.308534775    97.405272011      -3.4946846 %PBLI
        
        
    ];

filenames = [
        'BNON_2010-2013_postproc.mat'
        'BSIM_2010-2013_postproc.mat'
        'LEWK_2010-2013_postproc.mat'
        'BITI_2010-2013_postproc.mat'
        'BTHL_2010-2013_postproc.mat'
        'BSAT_2010-2013_postproc.mat'
        'LAIS_2010-2013_postproc.mat'
        'ABGS_2011-2013_postproc.mat'
        'BTET_201-0-1-3_postproc.mat'
        'SDKL_2012-2012_postproc.mat'
        'RNDG_2012-2012_postproc.mat'
        'PBKR_2012-2012_postproc.mat'
        'PBLI_2012-2012_postproc.mat'
    ];

gps_names = [
    'BNON'
    'BSIM'
    'LEWK'
    'BITI'
    'BTHL'
    'BSAT'
    'LAIS'
    'ABGS'
    'BTET'
    'SDKL'
    'RNDG'
    'PBKR'
    'PBLI'
    ];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process diurnal cycles for GPS stations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nf = size(filenames,1);

dcyc_mean_yr = zeros(24,Nf);
dcyc_stdm_yr = zeros(24,Nf);
norm_yr      = zeros(24,Nf);

for nf = 1:Nf
    
    load(filenames(nf,:))
    
    dcyc_mean_yr(:,nf) = dcyc_mean_tot(:);
    dcyc_stdm_yr(:,nf) = dcyc_stdm_tot(:);
    norm_yr(:,nf)      = norm_tot(:);
    
    clear dcyc_mean_tot
    clear dcyc_stdm_tot
    
end

dcyc_mean_yr(norm_yr < 100) = NaN;
dcyc_stdm_yr(norm_yr < 100) = NaN;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process diurnal cycle for COSMIC points %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ls = 0.5;
nf = 1;

pwv_cos = zeros(1,8);
std_cos = zeros(1,8);
num_cos = zeros(1,8);

for year = 2006:2012
    data_raw = importdata(['../data/cosmic_3hr/',num2str(year),'.txt']);
    
    Ndata = size(data_raw.data,1);
    ntmp  = 0;
    for nd = 1:Ndata
        
        if(...
                abs(data_raw.data(nd,5) - gps_chosen(nf,1)) < ls && ...
                abs(data_raw.data(nd,6) - gps_chosen(nf,2)) < ls ...
                )
            
            hr = floor(data_raw.data(nd,3)/3)+1;
            
            pwv_cos(hr) = pwv_cos(hr) + data_raw.data(nd,7);
            num_cos(hr) = num_cos(hr) + 1;
            
            ntmp = ntmp+1;
        end
        
    end
    
    disp([year ntmp])
    
end

pwv_cos = pwv_cos./num_cos;

for year = 2006:2012
    data_raw = importdata(['../data/cosmic_3hr/',num2str(year),'.txt']);
    
    Ndata = size(data_raw.data,1);
    ntmp  = 0;
    for nd = 1:Ndata
        
        if(...
                abs(data_raw.data(nd,5) - gps_chosen(nf,1)) < ls && ...
                abs(data_raw.data(nd,6) - gps_chosen(nf,2)) < ls ...
                )
            
            hr = floor(data_raw.data(nd,3)/3)+1;
            
            std_cos(hr) = std_cos(hr) + (data_raw.data(nd,7)-pwv_cos(hr))^2;
            
            ntmp = ntmp+1;
        end
        
    end
    
    disp([year ntmp])
    
end

std_cos = sqrt(std_cos./num_cos)./sqrt(num_cos);


pwv_gps = squeeze(dcyc_mean_yr(1:3:24,nf));
std_gps = squeeze(dcyc_stdm_yr(1:3:24,nf));



figure('Position', [1 1 800 800])
hold on
e1 = errorbar(0:3:21,pwv_cos,std_cos);
e1.LineWidth = 2;

e2 = errorbar(0:3:21,pwv_gps,std_gps);
e2.LineWidth = 2;
hold off

grid on
box on

xlim([0 21])

set(gca, 'xTick', 0:3:21)

xlabel('Time (hr)')
ylabel('PWV (mm)')

title(gps_names(nf,:))

legend('COSMIC', 'GPS')

