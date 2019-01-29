
nf    = 13;

Nt    = 1440;
Nz    = 44;

year  = 2011;



ti = 1;
if(nf == 1)
    ti = 481;
end

tf = 1440;
if(nf == 13)
    if(year == 2011)
        tf = 229;
        Nt = 229;
    end
    if(year == 2012)
        tf = 277;
        Nt = 277;
    end
end







load('wrf_domain.mat')

root  = '/n/regal/kuang_lab/gtorri/';
fn    = dir([root,'wrfout_d02_',num2str(year),'*']);
N_wrf = length(fn);


u10    = zeros(Nlat,Nlong,48);
v10    = zeros(Nlat,Nlong,48);
qv_sfc = zeros(Nlat,Nlong,48);
th_sfc = zeros(Nlat,Nlong,48);
p_sfc  = zeros(Nlat,Nlong,48);
sst    = zeros(Nlat,Nlong,48);

ncp    = netcdf.open([root,fn(nf).name],'NOWRITE');


tic
varid   = netcdf.inqVarID(ncp,'U10');
u_tmp   = netcdf.getVar(ncp, varid);

varid   = netcdf.inqVarID(ncp,'V10');
v_tmp   = netcdf.getVar(ncp, varid);

varid   = netcdf.inqVarID(ncp,'QVAPOR');
qv_tmp  = netcdf.getVar(ncp, varid, [0 0 0 0], [Nlat Nlong 1 Nt]);

varid   = netcdf.inqVarID(ncp,'T');
th_tmp  = netcdf.getVar(ncp, varid, [0 0 0 0], [Nlat Nlong 1 Nt]);

varid   = netcdf.inqVarID(ncp,'PSFC');
p_tmp   = netcdf.getVar(ncp, varid);

varid   = netcdf.inqVarID(ncp,'SST');
sst_tmp = netcdf.getVar(ncp, varid);
toc



for t = ti:tf
    
    ts = rem(t-1,48)+1;
    u10(:,:,ts)   = u10(:,:,ts)  + u_tmp(:,:,t);
    v10(:,:,ts)   = v10(:,:,ts)  + v_tmp(:,:,t);
    
    p_sfc(:,:,ts)   = p_sfc(:,:,ts)  + squeeze(p_tmp(:,:,t));
    qv_sfc(:,:,ts)  = qv_sfc(:,:,ts) + squeeze(qv_tmp(:,:,t));
    th_sfc(:,:,ts)  = th_sfc(:,:,ts) + squeeze(th_tmp(:,:,t));
    sst(:,:,ts) = sst(:,:,ts)  + sst_tmp(:,:,t);
    
    disp(t)
end

ndays = (tf-ti+1)/48;


u10   = u10/ndays;
v10   = v10/ndays;
p_sfc   = p_sfc/ndays;
qv_sfc  = qv_sfc/ndays;
th_sfc  = th_sfc/ndays;
sst = sst/ndays;

disp(nf)

save(['wrf_sfc_data_',num2str(nf),'.mat'], 'u10', 'v10', 'p_sfc', 'qv_sfc', 'th_sfc', 'sst')

