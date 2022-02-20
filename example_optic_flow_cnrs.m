
% Written by Sham Tlili and edited by Timothy Saunders
% Code for performing optic flow analysis


addpath('./useful_functions/')
% Defining outputs
screen_dpi          = get(0, 'ScreenPixelsPerInch');
screen_display_factor = 100;

% Input and output folders and files
folder_image        = './images/';
folder_save         = './results/flows/';
name                = 'kalman_filter.tif';
namebinar           = 'binar.tif';

% Threshold and size options
dgrid               = 10;
threshold           = 0.5;
median_threshold    = 2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%    KLT parameters     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.maxFeature    = 5000;     % feature number
param.winSize       = [24 24];  %window number
param.pyramid       = 2;        % pyramid size
param.maxIteration  = [200 200];% Iteration number
param.threshold     = [.1 .1 ]; % Thresholding
blurradius          = 1;        % Image blurring
 
%%%%%%%%%%%%%%%%%%%%%% General parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt                  = 1;        % time step  
istart              = 47;       % analysis start
iend                = 48;       % analysis end
length              = 1;        % arrow length 
scaleV              = 10;       % velocity scaling

%%%%%%%%%%%%%%%%%%%% start the code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define outputs
Xlag = NaN*ones(param.maxFeature,iend-istart+dt);
Ylag = NaN*ones(param.maxFeature,iend-istart+dt);
Ulag = NaN*ones(param.maxFeature,iend-istart+dt);
Vlag = NaN*ones(param.maxFeature,iend-istart+dt);
Time = NaN*ones(iend-istart+dt,1);

for i=istart:iend-dt
    %compute
    disp(['On number ', num2str(i), ' of ', num2str(iend)])  
    ind         = i-istart+1;
    Time(ind,1) = i;
    dtime = dt;                 %if needed changing time interval during time

    %create    
     if i== istart
        a=imread([folder_image name],i);
        [y,x]   = ndgrid(param.winSize:dgrid:size(a,1)-param.winSize,param.winSize:dgrid:size(a,2)-param.winSize);
        size_veul = size(x);  
        x_eul   = reshape(x,size_veul(1)*size_veul(2),1);
        y_eul   = reshape(y,size_veul(1)*size_veul(2),1);

        XEul    = NaN*ones(size_veul(1)*size_veul(2),iend-istart+dt);
        YEul    = NaN*ones(size_veul(1)*size_veul(2),iend-istart+dt);
        UEul    = NaN*ones(size_veul(1)*size_veul(2),iend-istart+dt);
        VEul    = NaN*ones(size_veul(1)*size_veul(2),iend-istart+dt);   
        index_points = NaN*ones(size_veul(1)*size_veul(2),iend-istart+dt);
        sa      = size(a);
        ysize   = sa(1);
        xsize   = sa(2);
     end

    a           = imread([folder_image name],i);  % open the analysis image
    
    acadre      = imread([folder_image namebinar],1);  % and the image mask
    acadre      = acadre/255;
    
    totalbinar  = acadre;
    
    a(find(totalbinar==0)) = NaN;

    for ss = 1:size(x_eul,1)
        box_values = totalbinar(y_eul(ss)-round(param.winSize(1)/2):y_eul(ss)+round(param.winSize(1)/2),x_eul(ss)-round(param.winSize(1)/2):x_eul(ss)+round(param.winSize(1)/2));
        box_values_mean(ss) = nanmean(nanmean(box_values));
        if box_values_mean(ss)>threshold
            index_points(ss,ind) = 1;
        else
            index_points(ss,ind) = 0;
        end
    end

    % Compute pyramid for first image (initially stored in second pyramid)
    pyr1 = makePyramid(a,param.pyramid,blurradius);

    % find local maxima            
    adouble = double(a);
    clear pts
    pts(:,2) = x_eul(find(index_points(:,ind)>0));
    pts(:,1) = y_eul(find(index_points(:,ind)>0));

    % open image for analysis
    b       = imread([folder_image name],i+dtime);
    bbinar  = imread([folder_image namebinar],1);
    bbinar  = 1-bbinar/255;
   
    totalbinar = acadre;
    
    b(find(totalbinar==0)) = NaN;
    pyr2    = makePyramid(b,param.pyramid,blurradius);
    [sp warn] = pyrLK(pyr1, pyr2, pts, ...
                param.winSize, param.maxIteration, param.threshold);
    warn    = sum(warn,2);
    indf    =~ warn;

    % First filter for NaN
    sizz    = size(pts(:,1));
    norm    = NaN*ones(sizz(1),1);
    
    sp      = sp/dtime;
    for ll = 1:sizz(1)
        norm(ll,1) = sqrt(sp(ll,1)* sp(ll,1)+sp(ll,2)* sp(ll,2));
    end

    normmean = nanmedian(norm);

     for ll = 1:sizz(1)
         disp(['On ', num2str(ll), ' of ', num2str(sizz(1))])

         XEul(find(index_points(:,ind)>0),ind) =  x_eul(find(index_points(:,ind)>0));
         YEul(find(index_points(:,ind)>0),ind) =  y_eul(find(index_points(:,ind)>0));
         UEul(find(index_points(:,ind)>0),ind) =  sp(:,2);
         VEul(find(index_points(:,ind)>0),ind) =  sp(:,1);

     end

     h12    = figure('PaperPositionMode','auto');
    imshow(a,[],'InitialMagnification',screen_display_factor,'Border','tight');
    hold all

    figure(h12)

    hold on
    hold all
    quiver(XEul(:,ind),YEul(:,ind),scaleV*UEul(:,ind),scaleV*VEul(:,ind),'color',[0 1 0],'linewidth',0.5,'Autoscale','off');
    set(gca,'ydir','reverse','Xtick',[],'Ytick',[])
    axis equal
    axis tight

    screen_dpi = get(0, 'ScreenPixelsPerInch');
    print(h12,'-dtiffn', sprintf('-r%d', 400*screen_dpi/screen_display_factor),  [folder_save,'V_Lag' num2str(i,'%03d') '.tif']);
    close(h12)

    resultats.XEul = XEul;
    resultats.YEul = YEul;
    resultats.UEul = UEul;
    resultats.VEul = VEul;

    resultats.x_eul= x_eul;
    resultats.y_eul= y_eul;

    resultats.Time = Time;
    resultats.dtime = dtime;
    resultats.dgrid = dgrid;

    clear pts
    clear sp
    close all 

    save([folder_save 'resultats.mat'],'resultats','-mat')

 end

