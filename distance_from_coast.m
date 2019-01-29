clear all

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultaxesfontsize', 19)

load('../data/wrf_domain.mat')


% Find the connected components of the figures and select the biggest
CC            = bwconncomp(lmask);
numPixels     = cellfun(@numel,CC.PixelIdxList);
[~,idx]       = max(numPixels);
Sumatra       = CC.PixelIdxList{idx};
[Sum_i,Sum_j] = ind2sub([Nlong,Nlat],Sumatra);


% Draw the island of sumatra
s_mask = zeros(Nlong,Nlat);
t_mask = zeros(Nlong,Nlat);
u_mask = zeros(Nlong,Nlat);
for n = 1:length(Sum_i)
    i = Sum_i(n);
    j = Sum_j(n);
    
    s_mask(i,j) = 1;
end


% Draw its contour
B       = bwboundaries(s_mask);
mat_bnd = B{1};
for n = 1:size(mat_bnd,1)
    i = mat_bnd(n,1);
    j = mat_bnd(n,2);
    
    t_mask(i,j) = 1;
end


% Select the bottom-right and top-left points
val_x = min(mat_bnd(:,2));
idi   = find(mat_bnd(:,2) == val_x);
val_y = max(mat_bnd(idi,1));
BR    = [val_y,val_x];

val_y = min(mat_bnd(:,1));
idj   = find(mat_bnd(:,1) == val_y);
val_x = max(mat_bnd(idj,2));
TL    = [val_y,val_x];


% Cut out the two extermes
for a1 = -1:2
    for a2 = -1:1
        i = TL(1);
        j = TL(2);
        t_mask(i+a1,j+a2) = 0;
        
        i = BR(1);
        j = BR(2);
        t_mask(i+a1,j+a2) = 0;
    end
end


% Find the smallest connected component of the boundary
DD            = bwconncomp(t_mask);
numPixels     = cellfun(@numel,DD.PixelIdxList);
[~,idx]       = min(numPixels);
coast         = DD.PixelIdxList{idx};
[cst_i,cst_j] = ind2sub([Nlong,Nlat],coast);
for n = 1:length(coast)
    i = cst_i(n);
    j = cst_j(n);
    u_mask(i,j) = 1;
end
dist_mat = bwdist(u_mask)*3;


% Make distances negative to the norht-east of the island
for j = 1:Nlat
    for i = 1:Nlong
        
        X = [i,j];
        Y = [cst_i,cst_j];
        ind = dsearchn(Y, X);
        tmp = Y(ind,:)-X;
        if(tmp(1) < 0 || tmp(2) < 0)
            dist_mat(i,j) = -dist_mat(i,j);
        end
        
    end
end

save('dist_mat.mat', 'dist_mat')

