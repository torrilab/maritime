
clear all
load('gps_data_2010-2012.mat')


Nf = size(gps_data,1);

para   = zeros(3,Nf);
dist   = zeros(1,Nf);
lats   = zeros(1,Nf);
longs  = zeros(1,Nf);
height = zeros(1,Nf);

for nf = 1:Nf
    
    
    name = gps_data(nf,1);
    coor = cell2mat(gps_data(nf,2));
    data = cell2mat(gps_data(nf,3));
    
    
    % Read the data in a format to feed the regression
    y    = data(:,8);
    x    = 2*pi*((data(:,1)-2010).*yeardays(data(:,1)) + data(:,2));
    
    yu   = max(y);
    yl   = min(y);
    yr   = (yu-yl);

    
    % Do the regression
    ym   = mean(y);
    fit  = @(b,x)  b(1).*(cos(x - b(2))) + b(3);
    fcn  = @(b) sum((fit(b,x) - y).^2);
    para(:,nf) = fminsearch(fcn, [yr;  0;  ym]);
    
    
    
end


% Normalize time lag
para(2,:) = para(2,:)/(2*pi)*24;


% Impose that amplitude be positive 
vec1 = para(1,:);
vec2 = para(2,:);

vec2(vec1 < 0) = vec2(vec1 < 0)+12;
vec1 = abs(vec1);

para(1,:) = vec1;
para(2,:) = vec2;





figure('Position', [1 1 800 800])
hold on
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
for nf = 1:Nf
    
    name = gps_data(nf,1);
    coor = cell2mat(gps_data(nf,2));
    
    scatter(coor(2),coor(1),100,para(2,nf),'fill');
    
end
hold off
box on
xlim([94 110])
ylim([-8 8])
xlabel('Longitude')
ylabel('Latitude')


