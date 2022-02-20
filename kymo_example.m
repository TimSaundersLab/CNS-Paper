addpath('/useful functions/');
velocity= importdata('./kymo_ex_data/resultats.mat');

dxpix   = 0.955;        % Pixel size in um
dttime  = 5;               % Time step in minutes
vscale  = (dxpix/dttime);

xlag    = velocity.XEul;
ylag    = velocity.YEul;
ulag    = vscale*velocity.UEul;
vlag    = vscale*velocity.VEul;

% xlag = velocity.XEul;
% ylag = velocity.YEul;
% ulag = velocity.UEul;
% vlag = velocity.VEul;

tstart  = 1;
tend    = 189;

%% Load mask
ImgName     = './kymo_ex_data/cns_binar.tif';

%%
xmin = nanmin(xlag);
xmax = nanmax(xlag);

dx_slide = 10;
dx_count = 50;

xtab = xmin:dx_slide:xmax;

for tt = tstart:tend
    disp(['On time ',num2str(tt),' of ',num2str(tend)])
   
    image_binar = imread(ImgName,tt);
    [y,x] = ndgrid(1:size(image_binar,1),1:size(image_binar,2));

    xlag_temp = xlag(:,tt);
    ylag_temp = ylag(:,tt);
    ulag_temp = ulag(:,tt);
    vlag_temp = vlag(:,tt);

    index_x = 1;

    for xx = xtab

        indices_temp = find((xlag_temp-(xx-dx_count/2)).*(xlag_temp-(xx+dx_count/2))<0);        % IS THIS THE PROBLEM
        xlag_temp_xx = xlag_temp(indices_temp);
        ylag_temp_xx = ylag_temp(indices_temp);
        ulag_temp_xx = ulag_temp(indices_temp);
        vlag_temp_xx = vlag_temp(indices_temp);

        xlag_temp_xx_cns = NaN*ones(size(xlag_temp_xx));
        ylag_temp_xx_cns = NaN*ones(size(ylag_temp_xx));
        ulag_temp_xx_cns = NaN*ones(size(ulag_temp_xx));
        vlag_temp_xx_cns = NaN*ones(size(vlag_temp_xx));

        for ll = 1:size(xlag_temp_xx,1)

            if image_binar(ylag_temp_xx(ll),xlag_temp_xx(ll)) == 255                      % IS THIS THE PROBLEM
                xlag_temp_xx_cns(ll) = xlag_temp_xx(ll);
                ylag_temp_xx_cns(ll) = ylag_temp_xx(ll);
                ulag_temp_xx_cns(ll) = ulag_temp_xx(ll);
                vlag_temp_xx_cns(ll) = vlag_temp_xx(ll);
            end  
        end

        xmean(tt,index_x) = nanmean(xlag_temp_xx_cns);
        ymean(tt,index_x) = nanmean(ylag_temp_xx_cns);
        umean(tt,index_x) = nanmean(ulag_temp_xx_cns);
        vmean(tt,index_x) = nanmean(vlag_temp_xx_cns);

        xstd(tt,index_x) = nanstd(xlag_temp_xx_cns);
        ystd(tt,index_x) = nanstd(ylag_temp_xx_cns);
        ustd(tt,index_x) = nanstd(ulag_temp_xx_cns);
        vstd(tt,index_x) = nanstd(vlag_temp_xx_cns);

        index_x= index_x+1;
    end

end

umean_smooth = tsmovavg_sham_gaussian_matrix(umean,[1 3 5 7 9 11 13 11 9 7 5 3 1],2);
umean_div = NaN*ones(size(umean_smooth));

for ll = 2:size(umean,2)-1
%     umean_div(:,ll) = (umean_smooth(:,ll+1)-umean_smooth(:,ll-1))/(2*dx_slide*dxpix);
    umean_div(:,ll) = (umean(:,ll+1)-umean(:,ll-1))/(2*dx_slide*dxpix);
end
    
%%

umean_smooth1   = umean_smooth + umean - umean;

%% Kymograph

time_involution = 1;
% dt  = 10;
min_velocity = -0.3;
max_velocity = 0.3;
kymo_t = NaN*ones(size(umean));
kymo_x = NaN*ones(size(umean));

for tt = 1:size(umean,1)
    kymo_t(tt,:) = (tt-time_involution)*dttime;
end

for xx = 1:size(umean,2)
    kymo_x(:,xx) = xtab(xx);
end
figure
% subplot(2,1,1)
hold on

g = surf(kymo_x,kymo_t,umean);
scatter3(100,0,1000,2,[0 0 0],'filled')
scatter3(100,600,1000,2,[0 0 0],'filled')
scatter3(500,0,1000,2,[0 0 0],'filled')
scatter3(500,600,1000,2,[0 0 0],'filled')
set(g, 'edgecolor','none')
% plot([kymo_x(1,20),kymo_x(1,20)],[1,600],'k-','LineWidth',4)
view(2)
set(gca,'FontSize',20)
h = colorbar;
colormap(fire)
caxis([min_velocity max_velocity])
hold off
