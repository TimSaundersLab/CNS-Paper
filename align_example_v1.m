
% Written by Sham Tlili and edited by Timothy Saunders
% Code for aligning data of CNS from SPIM data to detwich embryo
% Reads data in "images_original"
% tstart= first time point to be aligned
% tend = final time point to be aligned
% dbin = binning of data
% optimizer allows control of iterations, step lengths 

addpath('./Useful Functions/')
tic

folder_image    = './images_original/';     % Input data
folder_results  = './images_aligned/';      % Output data

tstart          = 1;                        % Start time
tend            = 2;                        % End time

dbin            = 8;                        % Binning size


for time = tstart+1:1:tend
    
    disp(['On time ',num2str(time), ' of ', num2str(tend)])
    % Load data
    name_movie1 = [folder_results 'image_' num2str(time-1) '.tif'];
    stack1      = double(readTIFstack(name_movie1));
    name_movie2 = [folder_image 'image_' num2str(time) '.tif'];
    stack2      = double(readTIFstack(name_movie2));

    [sy,sx,sz]  = size(stack1); 
    
    [y_pixels, x_pixels, z_pixels] = ndgrid(1:1:sy,1:1:sx,1:1:sz); 
    [y_bin, x_bin, z_bin]          = ndgrid(1:dbin:sy,1:dbin:sx,1:dbin:sz); 
    
    bined_stack1 = interp3(stack1,x_bin,y_bin,z_bin);
    bined_stack2 = interp3(stack2,x_bin,y_bin,z_bin);

    [optimizer, metric]         = imregconfig('monomodal');
    optimizer.MaximumIterations = 2000;
    %   optimizer.MinimumStepLength = 0.05;         % Can be adjusted if needed
    %   optimizer.MaximumStepLength = 0.1;          % Can be adjusted if needed

    tform       = imregtform(bined_stack2, bined_stack1, 'rigid', optimizer, metric);
    tform_final = tform;
    tform_final.T(4,1:3) = dbin*tform_final.T(4,1:3);

    slice_interp_turn    = imwarp(stack2,tform_final,'OutputView',imref3d(size(stack2)));
     
    name_align  = [folder_results 'image_' num2str(time) '.tif'];
    stack_turned = uint16(slice_interp_turn);
    writeTIFstack(stack_turned, name_align,  2^31)
   
end

toc


