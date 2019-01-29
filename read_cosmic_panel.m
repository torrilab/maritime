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


lat_1 = -4;
lat_2 =  4;

lon_1 =  95;
lon_2 = 103;

nt = 4;

figure('Position', [1 1 24000 1200])

for nt = 1:8
    
    subplot(2,4,nt)
    
    hold on
    geoshow('landareas.shp', 'EdgeColor', 'k', 'FaceColor', 'none');
    for year = 2006:2012
        data_raw = importdata(['../data/cosmic_3hr/',num2str(year),'.txt']);
        
        Ndata = size(data_raw.data,1);
        ntmp  = 0;
        for nd = 1:Ndata
            
            hr = floor(data_raw.data(nd,3)/3);
            lat_cos = data_raw.data(nd,5);
            lon_cos = data_raw.data(nd,6);
            pwv_cos = data_raw.data(nd,7);
            
            if(lat_cos > lat_1 && lat_cos < lat_2 ...
                    && lon_cos > lon_1 && lon_cos < lon_2 ...
                    && (hr == nt-1))
                
                scatter(lon_cos,lat_cos,150,pwv_cos, 'fill', 'markeredgecolor', 'k')
                
                ntmp = ntmp+1;
            end
            
        end
        
        disp([year ntmp])
        
    end
    
    
    for nf = 1:Nf
        scatter(gps_chosen(nf,2),gps_chosen(nf,1),300, dcyc_mean_yr(3*(nt-1)+1,nf) ,'d','fill', 'markeredgecolor', 'k')
    end
    
    box on
    
    caxis([30 60])
    colorbar
    
    xlim([95 103])
    ylim([-4 4])
    
    xlabel('Latitude')
    ylabel('Longitude')
    
    cb = colorbar;
    title(cb, '(mm)')
    
    title([num2str(3*nt),' hours'])
    
end

