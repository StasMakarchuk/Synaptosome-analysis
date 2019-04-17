%% Synapto-Analysis

% This script reads in a set of reconstructed single molecule localisation
% data sets from three different color channels and it automatically
% detects individual synaptosomes from each reconstructed field of view. 
% The script performs simple filtering on the data sets based on frame
% number, sigma value of the Gaussian fit, localisation precision, and
% intensity. A density filter is implemented to remove background
% localisations and only keep clustered detections. 
%
% The script outputs:
% - Reconstructed images for the three different color channels 
% - Overlay image of the three channels 
% - Binary image with potential detected synaptosomes in the red channel 
% - Text file with a summary of the thresholds established, and potential
%   number of synaptosomes detected. 
% - CSV file with centroids and weighted overlap calculations for each
%   potentially detected synaptosome. 

% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated: 10.04.19, Stas Makarchuk - extended for any number of
% channels


clear all
close all
clc

tic

%% File parameters

%directory = 'F:\synaptosomes\2018_10_10_Pedro_5thRound_EGTAK\output_reconstructions\Registered_data';
%directory = 'E:\Experiments\synaptosomes\Datasets_synaptosomes_20181206_4C_37C\37Cb\Data\thunderSTORM_phys\reconstructions\Registered_data';
directory = 'D:\PostDoc Cambridge\STORM\synapto_data_from_Janin\Sample\output_reconstructions\Registered_data';

% path to folder where outputfolder will be created (if doesn't already exist)
output_dir = fullfile('D:\PostDoc Cambridge\STORM\synapto_data_from_Janin\Sample\output_reconstructions\Registered_data',filesep);

repeat           = '4Cb';
condition        = 'phys';
%Specify here channels names that are used for file namings (could be one or more channels)
channel_token    = {'_RC', '_GC_reg'};
N_channel        = size(channel_token,2); 

% Single molecule reconstruction settings

pixelsize            = 117; % pixelsize in nm
magnification        = 10; % value of 10 gives 11.7 nm pixels in reconstruction (if pixelsize camera is 117 nm)
show                 = 1; % 1 to show extra intermediate results
show_d               = 0;
format               = 'thunderstorm'; % reconstruction software used (only thunderstorm)
format_for_filtering = 'thunderstorm';

% Filtering parameters

max_radius     = [0; 0]; % maximum radius of clusters in all chanels respectively (should be the same size as channel_token)
min_nr_locs    = [0; 0];  % minimum nr of locs within min_radius

P              = [0; 0]; % blobs with fewer pixels will be removed from mask


level          = [0; 0] % threshold (0-1) for mask (0 for Otsu's threshold)
RefCh          = 1;     %number of the channel which is areference channel for synaptosome analysis

sigma_kernel   = 15; % sigma for generating an image (if format = rapidstorm)
remove_border  = 0; % to remove the edges of images if there are edge artefacts
border_width   = 700; % width of border removed (in image coordinates)

filter          = 1;    % to filter the localisation files before further processing based on the following parameters
min_frames      = 0;  % throw away initial frames to avoid really high density data
min_sigma       = 40;   % maximum sigma value of the Gaussian fit
max_sigma       = 400;  % maximum sigma value of the Gaussian fit
min_intensity   = 1000; % min number of photons expected for an event
max_uncertainty = 40;   % max uncertainty expected for single molecules in focus
min_uncertainty = 5;    % max uncertainty expected for single molecules in focus 

essence = 0; % 1 to only write away files essential for further analysis
warning('off','images:initSize:adjustingMag');
warning('off','MATLAB:MKDIR:DirectoryExists');


%% Creating output folder and parameter file

% Create new output folder
path_output = fullfile(output_dir,['Results_' repeat ' ' condition]);
if exist(path_output, 'dir')
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    quest = '\fontsize{12}An output folder ''Results'' already exists. If you continue, data in this folder might be overwritten.';
    answer = questdlg(quest,'Message','Cancel','Continue',opts);
    if strcmp(answer,'Continue')
        mkdir(path_output);
    else
        return
    end
else
    mkdir(path_output);
end

% Write away workspace variables in a parameter file
save(fullfile(path_output,'parameters.mat'));
 max_radius = max_radius/magnification;


% Get list of locfiles of first channel
if ismac
    filelist = dir([directory strcat('/*',channel_token{1},'.csv')]);
elseif ispc 
    filelist = dir([directory strcat('\*',channel_token{1},'.csv')]);
end
    filelist
for i = 1:size(filelist,1)
    
    % Get filenames
    
    for nch = 1:N_channel
        if nch == 1
            filename{nch} = filelist(i).name;
            area_token = strsplit(filename{nch},'_'); area_token = area_token{1};
        else
            filename{nch,1} = strcat(area_token,channel_token{nch},'.csv');
        end
    end
 
    
    % Write some info to command window
    disp('##### Synaptosome analysis #####')
    disp(' ')
    disp(['Folder:             ' directory])
    for nch = 1:N_channel
        disp(['File for ' num2str(nch) ' channel:   ' filename{nch}])
    end
    
    
    % Open a summary file for writing
    summary_file = fopen(fullfile(path_output,strcat(area_token,'_summary.txt')),'wt');
    
    % Read in files
    for nch = 1:N_channel  locs(nch).channel = readLocFile(fullfile(directory,filename{nch}),format); end
    
    
    % Record number of localisations before filtering
    for nch = 1:N_channel  num_locs_prefilter(nch) = size(locs(nch).channel,1); end
    
    % Write away filtered localisations
    for nch = 1:N_channel  writeLocFile(locs(nch).channel,fullfile(path_output,strcat(area_token,'_locs_', channel_token{nch}, '_raw.csv')),format); end
  
    
    % Estimate Fov
%     Assumes there will be at least one localisation close to the edges
%     of the Fov in one of the three channels
    data_for_Fov=[];
    for nch =1:N_channel  data_for_Fov = [data_for_Fov; locs(nch).channel.x]; end
    for nch =1:N_channel  data_for_Fov = [data_for_Fov; locs(nch).channel.y]; end
    
    Fov = estimateFov(data_for_Fov, pixelsize);
    disp(['Estimated Fov:      ' num2str(Fov)])
    disp(' ')
    
    
    %% Optional removal of edge artefacts
    % Big high intensity blobs around the edges of images can throw off the
    % synaptosome detector, so a border with some specified width can be
    % cropped out before doing any further processing.
    
    % For removing edge artefacts
    % Edited by Stas Makarchuk 12.04.19
    if remove_border        
        
        for nch=1:N_channel
            sLocs = size (locs(nch).channel.x,1);
            nLocs = 1;
            for iL = 1:sLocs
                
                if locs(nch).channel.x(iL) > border_width & locs(nch).channel.y(iL) > border_width & locs(nch).channel.x(iL) < Fov*pixelsize - border_width & locs(nch).channel.y(iL) < Fov*pixelsize - border_width 
                    locs_filtered(nch).channel(nLocs,:) = locs(nch).channel(iL,:);
                    nLocs = nLocs+1;
                end
            end
        end
         
     else
        locs_filtered = locs;
     end
    
%     
%     %% Optional filtering
%     % The registered localisations filtered to keep only localisations
%     % detected after 500 frames (because these early frames are usually not
%     % sparse enough to assure good reconstruction), keep only localisations
%     % that have a sigma between 40 nm and 400 nm (to avoid grid artefacts
%     % and exclude localizations that are out of focus) and uncertainty
%     % smaller than 40 nm and larger than 5 nm (the expected resolution of the reconstructed
%     % dSTORM images). For the red channel the intensity threshold is set at
%     % 1000 photons, and for the green and blue channels at 500 photons. 
%     
    


    if filter
        % Filtering all channels ---------------------------------------------------
        
        for nch = 1:N_channel
            locs_filtered(nch).channel = locs_filtered(nch).channel(locs_filtered(nch).channel.frame       > min_frames ,:);
            locs_filtered(nch).channel = locs_filtered(nch).channel(locs_filtered(nch).channel.intensity       > min_intensity ,:);
            locs_filtered(nch).channel = locs_filtered(nch).channel(locs_filtered(nch).channel.uncertainty       > min_uncertainty ,:);
            locs_filtered(nch).channel = locs_filtered(nch).channel(locs_filtered(nch).channel.uncertainty       < max_uncertainty ,:);
            if ~strcmp(format_for_filtering,'rapidstorm')
                locs(nch).channel = locs_filtered(nch).channel(locs_filtered(nch).channel.sigma   > min_sigma  ,:);
                locs(nch).channel = locs_filtered(nch).channel(locs_filtered(nch).channel.sigma   < max_sigma ,:);
            end
        end
        
        
    end
%     
%     % Get number of localisations after filter

      for nch=1:N_channel  num_locs_postfilter(nch) = size(locs_filtered(nch).channel,1);  end

%     % calculate overall rejection rate (percentage) after localisation filters
      for nch=1:N_channel  overall_rejection(nch) = ((num_locs_prefilter(nch) - num_locs_postfilter(nch))/num_locs_prefilter(nch))*100;  end
      
      %write in txt files ratio of rejected points
      fprintf(summary_file,'Overall rejection rate after localisation filtering \n');
      for nch=1:N_channel fprintf(summary_file,[num2str(nch) '  channel: %2.f %%\n'],overall_rejection(nch));   end


%     %% Write away filtered localisations
      for nch=1:N_channel writeLocFile(locs_filtered(nch).channel,fullfile(path_output,strcat(area_token,['_locs_' num2str(nch) '_channel_filtered.csv'])),format);   end
        
%     %% Generate SMLM reconstructions (before density filtering)
%     
   
    % Red channel
      for nch = 1:N_channel 
            X(nch).channel = locs_filtered(nch).channel.x;
            Y(nch).channel = locs_filtered(nch).channel.y;
            if strcmp(format_for_filtering,'thunderstorm')
                if (locs_filtered(nch).channel.sigma == 0)
                    sigma(nch) = sigma_kernel;
                else
                    sigma(nch) = median(locs_filtered(nch).channel.sigma)/magnification;
                end 
            elseif strcmp(format_for_filtering,'rapidstorm')
                sigma(nch) = sigma_kernel;
            end
            
            intensities = locs_filtered(nch).channel.intensity;
            img(nch).channel = generateImage(X(nch).channel,Y(nch).channel,sigma(nch),intensities,Fov,pixelsize,magnification);
            filename{nch} = char(strcat(area_token,['_reconstruction_' num2str(nch) '_channel.tif']));
            if ~essence; imwrite(flip(flip(img(nch).channel,2),1), fullfile(path_output, filename{nch})); end
            if show; figure; imshow(flip(flip(img(nch).channel,2),1)); title(['Reconstruction of ' num2str(nch) ' channel']); end
      end

%     % Get two-colour reconstruction red-green (only for visualization)
      if N_channel>1
          for nch1 = 1:N_channel-1
              for nch2 = nch1+1:N_channel
                    merged_reconstruction = uint8(zeros(size(img(nch1).channel,1),size(img(nch2).channel,2),3));
                    merged_reconstruction(:,:,1) = uint8(img(nch1).channel);
                    merged_reconstruction(:,:,2) = uint8(img(nch2).channel);
                    filename{nch} = char(strcat(area_token,['_reconstruction_of_' num2str(nch1) '_and_' num2str(nch2) '_channel.tif']));
                    if ~essence; imwrite(flip(flip(merged_reconstruction,2),1), fullfile(path_output, filename{nch})); end
                    if show; figure; imshow(flip(flip(merged_reconstruction,2),1)); title(['Reconstruction ' num2str(nch1) ' and ' num2str(nch2) ' channel']); end
%     
              end
          end
      end

%     % Get all colour reconstruction (only for visualization)
      if N_channel>2
          merged_reconstruction = uint8(zeros(size(img(1).channel,1),size(img(2).channel,2),3));
          for nch = 1:N_channel  merged_reconstruction(:,:,nch) = uint8(img(nch).channel); end
          filename = char(strcat(area_token,['_reconstruction_of_all_channels.tif']));
          if ~essence; imwrite(flip(flip(merged_reconstruction,2),1), fullfile(path_output, filename)); end
          if show; figure; imshow(flip(flip(merged_reconstruction,2),1)); title(['Reconstruction_of_all_channels']); end
      end
      
    
%     
%     %% Density filtering
%     % To improve detection of synaptosomes in the red channel, some density
%     % filtering is applied (the images are cluttered). Only localisations
%     % that have a minimum number of neighbouring localisations within some
%     % specified search radius are kept. Considering the radius of a
%     % synaptosome is about 500 nm, this value is chosen for the search
%     % radius. By looking at a number of images and testing out different
%     % parameters, the minimum number of localisations within this radius
%     % was chosen to be 500. The variables containing the unfiltered
%     % localisations are cleared after filtering to save memory.
%     
    fprintf(summary_file,'Rejection Rate after density filtering \n');
    % Density filtering all channels -------------------------------------------
    for nch = 1:N_channel
        if max_radius(nch)~= 0 && min_nr_locs(nch) ~= 0
            % Perform density filtering
            [locs_density_filtered(nch).channel,indeces(nch).channel,~] = ...
                nearestNeighbourDensityFilter([locs_filtered(nch).channel.x locs_filtered(nch).channel.y], ...
                max_radius(nch),min_nr_locs(nch),show_d); disp(' ');
            if isempty(locs_density_filtered)
                disp(['All localizations in the ' num2str(nch) ' channel were filtered out during density filtering!']);
                continue
            else
                nr_locs_before(nch) = size(locs_filtered(nch).channel.x,1);
                nr_locs_after(nch)  = size(locs_density_filtered(nch).channel,1);
                fprintf(summary_file,[num2str(nch)  ' channel: %1.f %%\n'], ((nr_locs_before(nch)-nr_locs_after(nch))/nr_locs_before(nch))*100);
            end
            % Generate image from filtered localizations
            X(nch).channel = locs_density_filtered(nch).channel(:,1);
            Y(nch).channel = locs_density_filtered(nch).channel(:,2);
            if strcmp(format_for_filtering,'thunderstorm')
                sigma(nch) = median(locs_filtered(nch).channel.sigma)/magnification;
            elseif strcmp(format_for_filtering,'rapidstorm')
                sigma(nch) = sigma_kernel;
            end
            %sigma_RC = sigma_kernel;
            intensities = array2table([indeces(nch).channel locs_filtered(nch).channel.intensity],'VariableNames',{'id','intensity'});
            intensities_filtered = intensities(intensities.id == 1,:);
            intensities_filtered = intensities_filtered.intensity;

            img(nch).channel = generateImage(X(nch).channel,Y(nch).channel,sigma(nch),intensities_filtered,Fov,pixelsize,magnification);
            filename{nch} = char(strcat(area_token,['_reconstruction_' num2str(nch) '_channel_density_filtered.tif']));
            if ~essence; imwrite(flip(flip(img(nch).channel,2),1), fullfile(path_output, filename{nch})); end
            if show; figure; imshow(flip(flip(img(nch).channel,2),1)); title(['Reconstruction of ' num2str(nch) ' channel after density filtering']); end
        else
            locs_density_filtered(nch).channel = [locs_filtered(nch).channel.x locs_filtered(nch).channel.y];
            fprintf(summary_file,[num2str(nch) ' channel: 100 %%\n']);
        end
    end
    
        clear var locs % Clear variable containing unfiltered data to save memory
%
%     
%     %% Get binary mask for all three channels
%     % First an image is generated from the filtered localisations. The
%     % magnification should be chosen to match the resolution of the
%     % reconstruction (in our case the pixelsize of the camera is 117 nm and
%     % the mean resolution of the reconstructions 30 nm (as determined by
%     % FRC), so a magnification = 10 is used to get effective pixelsize of
%     % 11.7 nm which meets the Nyquist criterion). The generated image is
%     % thresholded using Otsu's threshold to get a binary mask. In the red
%     % channel, blobs in the mask that are smaller than some amount of
%     % pixels P are filled (because they probably are not synaptosomes).
%     
    fprintf(summary_file,'Threshold (Otsu''s method or user-specified):\n');
%     
    % Binarize image to get a mask in all channels
     for nch = 1:N_channel   
        if level == 0 % Get automatic threshold
            level_i = graythresh(img(nch).channel);
            disp(['Automatic threshold for ' num2str(nch) ' channel: ' num2str(level_i)]);
        else
            level_i = level(nch);
            disp(['Specified threshold for red ' num2str(nch) ' channel: ' num2str(level_i)]);
        end
        fprintf(summary_file,[ num2str(nch) ' channel: %f\n'],level_i);
        mask(nch).channel = im2bw(img(nch).channel,level_i);

        % Remove blobs that are too small
        mask(nch).channel = bwareaopen(mask(nch).channel,P(nch));
        if show
            figure;
            imshow(flip(flip(mask(nch).channel,2),1));
            title(sprintf([num2str(nch) ' channel after bwareaopen (P = %d)'], P(nch)));
        end
        if ~essence
            filename{nch} = strcat(area_token,'_detected_synaptosomes.tif');
            imwrite(flip(flip(mask(nch).channel,2),1), fullfile(path_output,filename{nch}));
        end
     end


%     
%     % Get all channels image mask ------------------------------------------------------
    merged_mask = zeros(size(mask(1).channel,1),size(mask(1).channel,2),3);
    
    for nch = 1:N_channel merged_mask(:,:,nch) = mask(nch).channel; end
    Filename = char(strcat(area_token,'_masks_all_channels.tif'));
    imwrite(flip(flip(merged_mask,2),1), fullfile(path_output, Filename));
    if show; figure; imshow(flip(flip(merged_mask,2),1)); title('Merged channels'); end
%     
%     
%     %% Calculate overlap between channels for all potential synaptosomes
%     if we have more than one channel
    if N_channel>1
    %     % Get labeled mask of reference for synaptosome channel 
         mask_labeled = bwlabel(mask(RefCh).channel);
    %     
    %     % Get number of potential synaptosomes detected
         numberOfClusters = max(max(mask_labeled));
         disp([num2str(numberOfClusters) ' potential synaptosomes detected.']); disp(' ');
         fprintf(summary_file,'%d potential synaptosomes detected.\n',numberOfClusters);

        if numberOfClusters > 0

            % Get centroids of the synaptosomes
            synaptosomeMeasurements = regionprops(mask_labeled,mask_labeled,'all');
            centroids = [synaptosomeMeasurements.Centroid];
            x_centroid = centroids(1:2:end-1);
            y_centroid = centroids(2:2:end);

            % Get area of the synaptosomes
            areas = [synaptosomeMeasurements.Area];
            % areas = areas*(magnification^2); % correct here for magnification factor, or do it later in analysis

            % Loop over detected synaptosomes and calculate their overlap with clusters
            % in the green and blue channel as a percentage, and a weighted overlap
            % using the image intensities
            
            
            overlap_perc = zeros(N_channel, numberOfClusters); 
            weighted_overlap = zeros(N_channel, numberOfClusters); 
            for j = 1:numberOfClusters

                disp(['Synaptosome ' num2str(j) '/' num2str(numberOfClusters)]);

                % Get mask for current synaptosome only
                maskCurrentSynaptosome = mask_labeled;
                maskCurrentSynaptosome(maskCurrentSynaptosome~=j) = 0;
                maskCurrentSynaptosome(maskCurrentSynaptosome==j) = 1;

                % Calculate percentage area overlap with green channel
                % we take all channels excep
                for nch = 1:N_channel
                    overlap_with_mask(nch).channel = and(logical(mask(nch).channel),logical(maskCurrentSynaptosome));
                    overlap = (sum(overlap_with_mask(nch).channel)/sum(logical(maskCurrentSynaptosome)))*100;
                    overlap_perc(nch,j) = overlap;
                    disp([num2str(overlap) ['% overlap with ' num2str(nch) ' channel.']]);
                end
               
                % Apply mask to all channels
                for nch = 1:N_channel
                    if nch==RefCh
                        masked_SR_img(nch).channel   = logical(maskCurrentSynaptosome).*double(img(nch).channel);
                    else
                        masked_SR_img(nch).channel = logical(mask(nch).channel).*double(img(nch).channel);
                    end
                end
                
                % Get weighted overlap for all channels except reference
                % channel
                for nch = 1:N_channel
                    if nch~=RefCh
                    X(nch).channel = getColocCoefficient(masked_SR_img(RefCh).channel,masked_SR_img(nch).channel);
                    weighted_overlap(nch,j) = X(nch).channel.mandersCoeff;
                    end
                end
                
                disp(' ');
            end

            % Summarize meaurements in a table
            synaptosome_ID = 1:numberOfClusters;
            results = [synaptosome_ID' x_centroid' y_centroid' areas'];
            %go through all non ref channels for overlap, weighted overlap
            %and names
            Names=[];
            for nch=1:N_channel 
                if nch~=RefCh    
                results = [results overlap_perc(nch,:).' weighted_overlap(nch,:).']; 
                Names = {Names, ['Overlap_with_' num2str(nch) '_channel'], ['Weighted_overlap_with_' num2str(nch) '_channel']};    
                end
            end
            Names = Names(2:end);
             
            results = array2table(results,'VariableNames', [{'ID','xCentroid','yCentroid', 'Area',},Names]);
            filename2 = char(strcat(area_token,'_results.csv'));
            writetable(results,fullfile(path_output, filename2));

            fclose(summary_file);
            %% Plot results
            
            % Overlap the reference channel with all non -reference
            for nch=1:N_channel
                if nch~=RefCh
                    mask2channels = zeros(size(mask(RefCh).channel,1),size(mask(RefCh).channel,2),3);
                    mask2channels(:,:,1) = mask(RefCh).channel;
                    mask2channels(:,:,2) = mask(nch).channel;
                    if show
                        figure;
                        imshow(flip(flip(mask2channels,2),1));
                        title([num2str(RefCh) ' channel and ' num2str(nch) ' channel'],'FontSize',15)
                    end
                    if ~essence
                        filename3 = char(strcat(area_token,['_overlap_' num2str(nch) '_channel.tif']));
                        imwrite(flip(flip(mask2channels,2),1), fullfile(path_output, filename3));
                    end
                end
            end
           
        else
            continue
        end
    end
end

toc
