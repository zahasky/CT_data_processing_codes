% Raw_CT_FIJI_format_v2.m
% Christopher Zahasky
% 1/15/2016, updated 4/19/2018
% This code loads stack of images in clincal CT scan, centers the image,
% coarsens the image and crops the edges outside the core assuming a
% cyclindrical core. Read options carefully to determine proper input. Note
% that paths to the java files '\mij.jar' and '\ij.jar' must be described
% for this code to work. These files are located in the folder called
% 'Fiji.app'. Also note that files for each scan/series must be in separate
% folders. Use script 'move_ct_files.m' to assist with this.

clear all
close all

%%%%% BEGIN INPUT %%%%%%%%%%%%%%%%%%%%
% Path must be formated so that \\ between folders for FIJI
path_to_files = 'H:\\Experiment Data\\Deformation_band_cores_Strathclyde\\BS21\\CT_data\\BS21_imb_2';
% Prefix to be added to folders, may leave blank
save_file_name = '';

% Folder prefix name. For example if the file names are 'scan_12' then the
% folder prefix is 'scan_'
folder_prefix = '';
% Crop to inlet and outlet slices. Open stack in FIJI and determine first
% slice of core and last slice of core, so that coreholder is cropped out
% April 2018 BS21
inlet_slice = 5;
outlet_slice =82;

% scan number (i.e. 1,2,3
scan_numbers = [1:3];

% Slice width of CT scans 
z_vox_size = 1.25; % [mm]
xy_vox_size = 160/512; % [mm]

% Crop data 1= yes, 0 = no
crop_data = 1;

% show aquisition time, 1=yes, 0=no
show_time = 1;

% If you are averageing multiple scans then set this to 1, things are
% changing between scans then MAKE SURE THIS IS 0.
% average scans 1= yes, 0 = no
average_scans = 0;

% coarsen scans 1= yes, 0 = no
coarse_scans = 1;
% coarsen factor (i.e how many voxels are averaged in each direction). If
% you want voxel size of 1.25 mm then it should be cfx = 4 and cfz = 1. If
% you want voxel sizes of 2.5 mm then it should be cfx = 8 and cfz = 2.
% x direction
cfx = 8;
% y direction
cfy = cfx ;
% z direction
cfz = 2 ;
% width of core (in voxels). This determines the width of the core crop
% dimensions. Note that this number MUST BE DIVISBLE BY 'cfx'.
width_and_height = 144;

% Which slice to plot and compare
plot_slice = 6;

% FIJI input
% Setup paths to java files, this needs to be changed for whatever computer
% you are running this on
javaaddpath 'M:\Fiji.app\mij.jar'
javaaddpath 'M:\Fiji.app\ij.jar'

%%%%%%% END INPUT

cd([path_to_files, '\\', folder_prefix, num2str(scan_numbers(1))])

% files in folder
file_list = dir();
first_slice = file_list(3).name;

number_of_files = length(scan_numbers);
number_of_slices = outlet_slice-inlet_slice+1;
length_along_core = ([1:number_of_slices]*z_vox_size - z_vox_size/2)./10;

% Open FIJI
MIJ.start
figure

% Loop through images
for i = scan_numbers
    stack_folder = [folder_prefix, num2str(scan_numbers(i))];
    image_load_command = ['open=[', path_to_files, '\\', stack_folder,'\\', first_slice, '] sort'];
    
    % Load images as stacks
    MIJ.run('Image Sequence...', image_load_command);
    MIJ.run('32-bit');
    if show_time == 1
        MIJ.run('Show Info...')
    end
    
    % Begin core centering algorithm here and adjust crop box
    MIJ.run('Duplicate...', 'duplicate range=40-40')
    
    % apply threshold (this may need to be adjusted depending on the scans)
    MIJ.setThreshold(1100, 3500);
    MIJ.run('Convert to Mask');
    MIJ.run('Set Measurements...', 'mean standard center limit add redirect=None decimal=3');
    MIJ.run('Analyze Particles...', '  show=Masks display clear');
    x = MIJ.getResultsTable;
    % center of mass of core
    xm = (x(3,3)./xy_vox_size);
    ym = (x(3,4)./xy_vox_size);
    xcoord = round(xm - width_and_height/2);
    ycoord = round(ym - width_and_height/2);
    
    
   image_box = ['width=', num2str(width_and_height),' height=', ...
       num2str(width_and_height),' x=', num2str(xcoord), ...
       ' y=', num2str(ycoord), ' ,slice=10'];
    
    MIJ.selectWindow(num2str(stack_folder));

    MIJ.run('Specify...', image_box);
    
    MIJ.run('Crop');
    
    % If you want to interpolate data then uncomment this
    % MIJ.run('Size...', 'width=20 height=20 depth=96 constrain average interpolation=Bilinear');
    
    % convert image data to matlab matrix named 'CI'
    CI = MIJ.getCurrentImage();
    MIJ.closeAllWindows();
    
    % crop off ends of core
    CI = CI(:,:,inlet_slice:outlet_slice);
    
    % If you want to coarsen the scans (which you should), then do some
    % linear algebra
    if coarse_scans == 1
        crop_size = size(CI);
        
        % average original rows down into original/cfx
        v=repmat({ones(cfx,1)/cfx},1, floor(crop_size(1)/cfx));
        A = blkdiag(v{:});
        
        CI_coarse = zeros(crop_size(1)/cfx, crop_size(2)/cfy, crop_size(3)/cfz);
        for k = 1:crop_size(3)/cfz
            sum_slices = zeros(floor(crop_size(1)/cfx),floor(crop_size(2)/cfy));
            for j = 1:cfz
                slice = squeeze(squeeze(CI(:,:,k*cfz+j-cfz)));
                smaller_slice=A'*slice*A;
                sum_slices = sum_slices+smaller_slice;
            end
            CI_coarse(:,:,k) = sum_slices./cfz;
        end
        
        if crop_data == 1
            % crop outside of core and replace values with nans
            crop_dim = size(CI_coarse);
            r = (crop_dim(1)-1)/2;
            %         r = (crop_dim(1)+1)/2;
            cx = r+1;
            cy = r+1;
            dia = crop_dim(1);
            for ii=1:crop_dim(1)
                yp = round(cy + sqrt(r^2 - (ii-cx)^2));
                ym = round(cy - sqrt(r^2 - (ii-cx)^2));
                
                if yp <= dia
                    CI_coarse(ii,yp:end,:) = nan;
                end
                
                if ym >= 0
                    CI_coarse(ii,1:ym,:) = nan;
                end
            end
            cd(path_to_files)
            save([stack_folder,'_coarse.mat'], 'CI_coarse')
        elseif crop_data==0
            cd(path_to_files)
            save([stack_folder,'_coarse_uncropped.mat'], 'CI_coarse')
        end
        
    elseif coarse_scans == 0
        if crop_data == 1
            crop_dim = size(CI);
            % if a crop template is given use that
            if exist('crop_temp','var')
                for sl = 1:crop_dim(3)
                    cci = CI(:,:,sl);
                    cci(isnan(crop_temp)) = nan;
                    CI(:,:,sl) = cci;
                end
            else
            % crop outside of core and replace values with nans
            
            r = (crop_dim(1)-1)/2;
            %         r = (crop_dim(1)+1)/2;
            cx = r+1;
            cy = r+1;
            dia = crop_dim(1);
            for ii=1:crop_dim(1)
                yp = round(cy + sqrt(r^2 - (ii-cx)^2));
                ym = round(cy - sqrt(r^2 - (ii-cx)^2));
                
                if yp <= dia
                    CI(ii,yp:end,:) = nan;
                end
                
                if ym >= 0
                    CI(ii,1:ym,:) = nan;
                end
            end
            end
            cd(path_to_files)
%             save([stack_folder,'_cropped.mat'], 'CI')
        elseif crop_data == 0
            cd(path_to_files)
            save([stack_folder,'_raw.mat'], 'CI')
        end
    end
    
    % Plot slice to make sure it's right
    subplot(1,2,1)
    if coarse_scans == 1
        slice_exam = squeeze(squeeze(CI_coarse(:,:,1)));
    else
        slice_exam = squeeze(squeeze(CI(:,:,1)));
    end
    imagesc(slice_exam);
    axis equal
    axis tight
    shading flat
    colorbar
    title('coarsened front', 'fontsize', 14)
    subplot(1,2,2)
    if coarse_scans == 1
        slice_exam = squeeze(squeeze(CI_coarse(:,:,end)));
    else
        slice_exam = squeeze(squeeze(CI(:,:,end)));
    end
    imagesc(slice_exam);
    axis equal
    axis tight
    shading flat
    colorbar
    title('coarsened end', 'fontsize', 14)
    drawnow
    % If calculating an average of all of the images
    if average_scans == 1
        if i==1
            if coarse_scans == 1
                AVG = CI_coarse;
            else
                AVG = CI;
            end
        else
            if coarse_scans == 1
                AVG = AVG + CI_coarse;
            else
                AVG = AVG + CI;
            end
        end
    else
        SAVG = squeeze(squeeze(nanmean(nanmean(CI))));
        
        % plot(length_along_core, SAVG, 'k', 'linewidth', 2)
        % drawnow
        cd(path_to_files)
        save([save_file_name, num2str(scan_numbers(i)), '_sa_profile.mat'], 'SAVG')
    end
end

if average_scans == 1
    cd(path_to_files)
    % Average images
    AVG = AVG./number_of_files;
    
    % Plot slice to make sure it's right
    subplot(1,2,1)
    slice_exam = squeeze(squeeze(AVG(:,:,1)));
    imagesc(slice_exam);
    axis equal
    axis tight
    shading flat
    colorbar
    title('coarsened front', 'fontsize', 14)
    subplot(1,2,2)
    slice_exam = squeeze(squeeze(AVG(:,:,end)));
    imagesc(slice_exam);
    axis equal
    axis tight
    shading flat
    colorbar
    title('coarsened end', 'fontsize', 14)
    drawnow
    
    % Take slice averages
    SAVG = squeeze(squeeze(nanmean(nanmean(AVG))));
    save([save_file_name, '_sa_profile.mat'], 'SAVG')
    % Save as .mat file
    if coarse_scans == 1
        CI_coarse = AVG;
        save([save_file_name, '_coarse_average.mat'], 'CI_coarse')
    else
        CI = AVG;
        save([save_file_name, '_average.mat'], 'CI')
    end
end

% Close MIJ
MIJ.exit
