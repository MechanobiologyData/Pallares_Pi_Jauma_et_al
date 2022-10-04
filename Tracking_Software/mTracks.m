%% Software to segment phase constrast images
% Data are microscopy images stored in a stack named 'phas.tif'
% Images are located in 'dirBase'\StitchingPhase


clear all;
clc;

% Select folder were data is stored
dirBase = 'your_data_folder\';

% Set up folder where data will be stored
if ~exist([dirBase,filesep,'TrackingFiles']), mkdir([dirBase,filesep,'TrackingFiles']); end
if ~exist([dirBase,filesep,'SegmentationImages']), mkdir([dirBase,filesep,'SegmentationImages']); end

% Maximum number of frames in the ima stack
timeMax = 98;


%% Segmenting images into Black and White
% Results in "SegmentationImages" folder
% Adjust parameters for different images

% Parameters for filters
pixGAUSS = 9; 
pixMED = 1; 
pixSTD   = 4; 
pixGAUSS2 = 1;  
pixDIL   = 4;
pixEROD  = 18;

% Destination folders
dirInit = [dirBase, filesep, 'StitchingPhase'];
dirDest = [dirBase,filesep,'Tracking'];
    
if ~exist(dirDest), mkdir(dirDest); end

% Loop over time
for k=1:timeMax
    disp([' Time: ', num2str(k)]);

    % Open each frame of the image
    im = imread([dirInit, filesep, 'phas.tif'],k);

    % Normalize image
    if k==1
        meanIm0 = mean(double(im(:)));
    end
    meanIm = mean(double(im(:)));
    imNorm = (im/meanIm)*meanIm0;
    clear meanIm imCropped

    % 1st Gaussian Filter
    im2 = medfilt2(imNorm,[pixMED pixMED]);
    im2 = imgaussfilt(im2,pixGAUSS);

    % STD filtering
    SE = strel('disk',pixSTD);
    im2 = stdfilt(im2,SE.Neighborhood);

    % Threshold
    im2h = im2(:);
    im2h(im(:)==0) = [];
    [hist_img, bins] = imhist(im2h/max(im2(:)), 2^16);
    clear im2h

    % Triangle method for thresholding
    best_idx=triangle_th(hist_img,2^16);
    bw = imfill(im2bw(im2/max(im2(:)), best_idx), 8,'holes');

    % Erode
    bw2 = imerode(bw,strel('disk',pixEROD));
    %imagesc(im2), axis equal, colormap(gray);

    % Save image
    imSegPath = [dirBase,filesep,'SegmentationImages', filesep];
    if ~exist(imSegPath), mkdir(imSegPath); end
    imwrite(bw2,[imSegPath, filesep 'bw_', num2str(k),'.tif'] ,'tif', 'compression', 'none');
    clear im2 bw2

end

%% Tracking points

% Define the arrays to store data
points = {};
area = {};

% Define destination folders
dirInit = [dirBase, filesep, 'StitchingPhase'];
dirDest = [dirBase,filesep,'Tracking'];

% Loop to label objects within the image
for k=1:timeMax
    disp([' Time: ', num2str(k), ' ... Segmentation']);

    % Load image
    imSegPath = [dirBase,filesep,'SegmentationImages'];
    bw2 = imread([imSegPath, filesep 'bw_', num2str(k),'.tif']);

    % Remove objects touching the border
    bw2 = imclearborder(bw2);

    % Filter for small particles
    bwL = bwlabel(bw2);
    props = regionprops('table',bwL, 'Area', 'Centroid', 'Eccentricity');
    idx = find([props.Area] > median(props.Area) & [props.Eccentricity] < 0.95);
    bw3 = ismember(bwL, idx);

    % Store centroid and area of clusters
    points{k} = props.Centroid(idx,:);
    area{k} = props.Area(idx);

    clear bwL idx props
    clear im0 bw bw2 bw3 imC best_idx hist_img
end


points = points';
area   = area';

segmentation.pixGAUSS = pixGAUSS;
segmentation.pixSTD   = pixSTD;
segmentation.pixDIL   = pixDIL;
segmentation.pixEROD  = pixEROD;

% Save the data
save([dirDest, filesep, 'Points.mat'], 'points', 'area', 'segmentation');
save([dirBase,filesep,'TrackingFiles', filesep, 'Points','.mat'], 'points', 'area', 'segmentation');



%% Linking tracks
% So far we have individual points in each image. It's time to link those
% points to build individual trajectories

max_linking_distance = 200;
max_gap_closing = 1;
debug = false;

disp([' ... Linking']);


dirInit = [dirBase, filesep, 'StitchingPhase'];
dirDest = [dirBase,filesep,'Tracking'];

if ~exist(dirDest), mkdir(dirDest); end

load([dirBase,filesep,'TrackingFiles', filesep, 'Points','.mat']);

[ tracks, adjacency_tracks, A] = simpletracker(points,...
    'Method', 'NearestNeighbor',...%% 'Hungarian', ... %'NearestNeighbor',...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);


% Build TRACKS
clear tm
tm{numel(tracks)} = [];
for i = 1:numel(tracks)

    indices = tracks{i};
    values = [];

    for t=1:timeMax
        if ~isnan(indices(t))
            values = [values; t, points{t}(indices(t),:), area{t}(indices(t),:)];
        end
    end
    tm{i} = values;
    clear values t indices
end

% Store tracks
save([dirDest,filesep, 'Tracks.mat'], 'tm');
save([dirBase,filesep,'TrackingFiles',filesep, 'Tracks.mat'], 'tm');
clear points
    


%% Plot tracks to discard those poorly tracked

% Load tracks
load([dirBase,filesep,'TrackingFiles', filesep, 'Tracks.mat']);
colors = hsv(size(tm,2));
    
% Loop to create images of segmented clusters and their trajectories
for t=1:timeMax
    disp([' Time: ', num2str(t), ' ... Ploting for selection']);

    % Load phase image 
    im = imread([dirInit, filesep, 'phas.tif' ],t);

    % Load BW image
    imSegPath = [dirBase,filesep,'SegmentationImages'];
    bw = imread([imSegPath, filesep 'bw_', num2str(t),'.tif']);

    % Polygon around boundaries
    boundaries = bwboundaries(bw);

    clear bw imSegPath

    h= figure('visible','off');
    imagesc(im), axis equal, colormap(gray);
    hold on
    for j=1:size(boundaries,1)
        b = boundaries{j};
        plot(b(:,2),b(:,1),'LineWidth',1, 'Color', [232, 212, 255]/255);
        hold on;
    end
    hold on
    for i_track = 1:size(tm,2)

        if ~isempty(tm{i_track})
            time = tm{i_track}(:,1);
            ar = tm{i_track}(:,4);
            if length(time) > 0 & length(time(time==t))==1 & std(ar)/mean(ar) < 0.6 % Vatiation of area within the cluster

                valt = find(time==t);
                plot(tm{i_track}(valt,2),tm{i_track}(valt,3),'*', 'Color', colors(i_track, :))
                txt = num2str(i_track);
                text(tm{i_track}(valt,2),tm{i_track}(valt,3),txt,'Color', [0.98824, 0.9098,0.011765], 'FontSize',12);
                plot(tm{i_track}(1:valt,2), tm{i_track}(1:valt,3), 'Color', colors(i_track, :))

                clear valt
            end
        end
    end
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    axis off;
    export_fig(h, [dirDest,filesep,'plot_',num2str(t),'.png'])

    close all;
    end
    
