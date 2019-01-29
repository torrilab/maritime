



vec_names = zeros(0,4);
name      = char(1,4);


Nn        = 0;
gps_coo   = zeros(1,3);
gps_data  = cell(0,3);


for yr = 2008:2013
    
    
    % Read the coordinates of all the stations
    fn_pos   = [num2str(yr),'_Suominet_lat_lon_elevation.txt'];
    gps_pos  = importdata(['../data/gps/',fn_pos]);
    
    
    % Read the data for each station
    gps_dir  = dir(['../data/gps/',num2str(yr)]);
    Nf = length(gps_dir);
    
    
    for nf = 1:Nf
        
        N_data = 0;
        gps_1  = zeros(0,8);
        
        gps_name = gps_dir(nf).name;
        if(length(gps_name) > 2 && strcmp(gps_name(end-2:end),'dat')==1)
            
            namf = gps_name;
            name = gps_name(1:4);
            
            % Find position of the stations
            for np = 1:length(gps_pos.textdata)
                name_tmp = gps_pos.textdata{np};
                name_tmp = name_tmp(~isspace(name_tmp));
                
                if(strcmpi(name,name_tmp) == 1)
                    gps_coo(1) = gps_pos.data(np,1);
                    gps_coo(2) = gps_pos.data(np,2);
                    gps_coo(3) = gps_pos.data(np,3);
                    disp(name_tmp)
                end
            end
            
            % Import data
            gps_raw = importdata(['../data/gps/',num2str(yr),'/',namf]);
            pwv_raw = squeeze(gps_raw(:,9));
            id_r    = find(pwv_raw > 0);
            N_tmp   = length(id_r);
            
            if(N_tmp > 0)
                gps_tmp = zeros(N_tmp,8);
                for m = 1:N_tmp
                    n = id_r(m);
                    
                    gps_tmp(m,1) = gps_raw(n,1); % Year
                    gps_tmp(m,2) = gps_raw(n,2); % Time of year
                    gps_tmp(m,3) = gps_raw(n,3); % Month
                    gps_tmp(m,4) = gps_raw(n,4); % Day
                    gps_tmp(m,5) = gps_raw(n,5); % Hour
                    gps_tmp(m,6) = gps_raw(n,7); % Pressure
                    gps_tmp(m,7) = gps_raw(n,8); % Temperature
                    gps_tmp(m,8) = gps_raw(n,9); % PWV
                    
                end
                
                flag = 0;
                for nn = 1:Nn
                    if(size(vec_names,1) > 0 && contains(name,vec_names(nn,1:4)) == 1)
                        flag = 1;
                        break
                    end
                end
                
                if(flag == 1)
                    gps   = cell2mat(gps_data(nn,3));
                    gps_1 = cat(1,gps_tmp,gps);
                    gps_data(nn,1) = {name};
                    gps_data(nn,2) = {gps_coo};
                    gps_data(nn,3) = {gps_1};
                else
                    vec_names = cat(1,name,vec_names);
                    Nn    = Nn+1;
                    gps_data = cat(1,{name,gps_coo,gps_tmp},gps_data);
                    disp('nonexist')
                end
            end
            
            
            disp([name num2str(yr)])
        end
        
    end
    
end

save('gps_data_2008-2013.mat', 'gps_data');


vec   = gps_data(:,1);
[A,B] = sort(vec);
tmp   = gps_data(B,:);
gps_data = tmp;
save('gps_data_2008-2013_alph.mat', 'gps_data');
