%% Image Analysis
% location = input('Supply (in single quotes) the filepath to image folder');
% img = imread(location);
close all
location = '~/Desktop/College/Research/PayseurLab/bluemale.tif'; % DELETE
img = imread(location);
% % Creates list of all .tif files in the directory provided
% files = dir(strcat(location,'/*.tif')); %finds all matching files
% file_names = {files.name};
% file_names = file_names'; % transpose so correct dimensions
% file_path = cell(size(file_names,1), 1); % paste filepaths onto names
% file_path(:) = {location};
% file_paths = strcat(file_path,file_names);

flag = false; % DELETE - to mark when a bad image has been given
% in the case of a bad image, analysis will continue in spite of that by
% disregarding the previous expected values (e.g. centromere #)

% Determine if the slide is male or female from filepath
[~,name,~] = fileparts(location);
isFemale = regexp(name,'female');
isMale = isempty(isFemale);


% Split Channels -
%https://www.mathworks.com/matlabcentral/answers/91036-split-an-color-image-to-its-3-rgb-channels
red = img(:,:,1); % Red channel (BW)
green = img(:,:,2); % Green channel (BW)
blue = img(:,:,3); % Blue channel(BW)
a = zeros(size(img, 1), size(img, 2));
just_red = cat(3, red, a, a);
just_green = cat(3, a, green, a);
just_blue = cat(3, a, a, blue);

images = [img, just_red; just_green, just_blue];
figure, montage(images,'Size', [1 1]), title('Raw Color Channels');


%% Detection of Red

% enhance the red of the image - counter the gradual decrease in red
% intensity from left to right
stretch = decorrstretch(red,'tol',0.02);
darker = imadjust(stretch, stretchlim(stretch),[0.02 0.99]);

bw = imbinarize(darker, graythresh(darker)); % binarizes with best threshold
bwSharp = bwareaopen(bw, 50); % reduces background noise of binarized image
bwSmooth = imgaussfilt(double(bwSharp), 1); % makes the image smooth so that the fill method doesn't get confused
bwR = imbinarize(bwSmooth, graythresh(darker)); % binarizes again

%EDGE APPROACH - https://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
[~, threshold] = edge(bwR, 'sobel');
numFound = 0;
dilationFactor = 90;
fudgeFactor = .5;
dilationOn = 3; %little difference btwn 2 and 3, whereas 1 == off

BWs = edge(bwR,'sobel', threshold * fudgeFactor); % detect edge

% adjust dilation of lines from edge detection to connect loose fragments
while numFound ~= 20
    seBegin = strel('line', dilationOn, dilationFactor); % may combine very close chromosomes
    seEnd = strel('line', dilationOn, 0); % but also will combine fragmented ones
    BWsdil = imdilate(BWs, [seBegin seEnd]); % dilates edges
    final = imfill(BWsdil, 'holes'); % fill in area within edges
    
    redFound = bwconncomp(final, 8); % attempts to count number of objects
    numFound = redFound.NumObjects;
    
    % attempt to adjust
    if numFound < 20 && dilationFactor > 0 % too much dilation
        dilationFactor = dilationFactor - 5;
    elseif numFound > 20 && dilationFactor < 120 % too little dilation
        dilationFactor = dilationFactor + 5;
    elseif numFound < 20 && dilationFactor == 0 && dilationOn ~= 1 % turns off dilation completly
        dilationOn = 1;
    elseif numFound > 20 && dilationFactor == 120 && dilationOn < 5
        dilationFactor = 50;
        dilationOn = dilationOn + 1;
    elseif (numFound < 20 && dilationOn == 1) || dilationFactor >= 120
        % Dilation turned off, still can't do it... exit
        %         figure, imshowpair(img, label2rgb(labelmatrix(CC)), 'montage')
        %         title(strcat(['Failed, Chromosomes identified = ', num2str(numFound)]))
        warning('Bad image')
        beep
        flag = true;
        break
    end
    %     fprintf('NumFound: %i; Dilation Factor: %d; Dilation: %i\n',...
    %         numFound, dilationFactor, dilationOn) %DELETE
end

figure, imshowpair(img, label2rgb(labelmatrix(redFound)), 'montage')
title(strcat(['Binarized Red Channel, Chromosomes identified = ', num2str(numFound)]))

%% Detection of blue
% females have (19 autosomal + 1 sex) centromeres == 20
% males have (19 autosomal + 2 WEAK) centromeres == 19

if ~flag || redFound.NumObjects < 20 % second condition triggered by lots of overlaps
    numOfCentromeres = 20 - isMale; % searches for 20 if F, 19 if M
else % if bad image flag has been set, assume 1 chromosome = 1 centromere
    numOfCentromeres = redFound.NumObjects;
end

%defaults
numFound = 0;
threshold = graythresh(blue);
bckgrndReduct = 15; % to be adjusted
adjustCount = 0;
adjustLimit = 2000; % to catch infinite loops -> normal is < 100

while numFound ~= numOfCentromeres
    % Prepare Image for counting
    bw = imbinarize(blue, threshold); % binarizes with best threshold
    
    %Special Case: Auto-thresholding was very bad - arbitrarily make 0.7
    percentOfImg = sum(sum(bw))/ (size(img, 1) * size(img, 2));
    if(adjustCount == 0 && percentOfImg > 0.2)
        warning('special case thresholding: blue')
        threshold = 0.7;
        continue
    end
    
    bwBRaw = bwareaopen(bw, bckgrndReduct); % reduces background noise of binarized image
    bwB = bwareaopen(bwBRaw & bwR, bckgrndReduct); % overlap the red and blue channels to ensure
    
    %count number of objects found
    blueFound = bwconncomp(bwB, 8); % attempts to count number of objects
    numFound = blueFound.NumObjects;
    
    %adjust
    if numFound < numOfCentromeres && threshold > 0.005 % threshold is too high... not enough
        threshold = threshold - 0.005;
        adjustCount = adjustCount + 1;
    elseif numFound > numOfCentromeres && threshold < 0.995% threshold is too low, too many
        threshold = threshold + 0.005;
        adjustCount = adjustCount + 1;
    elseif threshold >= 0.995 %  max threshold used - adjust bckgrnd reduct.
        bckgrndReduct = bckgrndReduct + 1;
        threshold = graythresh(blue); %reset
    end
    
    if bckgrndReduct > 30
        warning('Repeated Exceptional Behavior - Continuing with best image')
        break
    end
    
    if adjustCount == adjustLimit
        warning('Exiting infinite adjustment loop - Continuing with best image')
        break
    end
end

figure, imshowpair(img, label2rgb(labelmatrix(blueFound)), 'montage')
title(strcat(['Binarized Blue Channel, Centromeres identified = ', num2str(numFound),...
    ', threshold = ', num2str(threshold)]))

%% Detection of green
% Legitimate green points of interest are very small - easily confused with
% background noise - must try different approach:
%  eliminate all areas of green that don't overlap red
bw = imbinarize(green, graythresh(green)); % binarizes with best threshold

% Special case: The foci are so faint that a bad threshold is selected
percent_SC_covered = sum(sum(bw & bwR))/ sum(sum(bwR));
if(percent_SC_covered > 0.3)
    warning('special case thresholding: green')
    bw = imbinarize(green, 0.175); %arbitrary
end

bw = bwareaopen(bw, 4); % reduces background noise of binarized image

% find foci
[~, threshold] = edge(bw, 'sobel');
bw = edge(bw,'sobel', threshold * fudgeFactor);

% enlarge
seBegin = strel('line', 2, 45);
seEnd = strel('line', 2, 0);
lines = imdilate(bw, [seBegin seEnd]);
bw = imfill(lines, 'holes');

% prune away foci that don't overlap red
bwG = bwareaopen((bwR & bw), 8);

greenFound = bwconncomp(bwG, 8); % attempts to count number of objects
numFound = greenFound.NumObjects;

figure, imshowpair(img, bwG, 'montage'),
title(strcat(['Binarized Green Channel, Foci identified = ', num2str(numFound)]))

%% Review the New Composite
% Simplify the image by only green and blue areas that overlap red
overlay = cat(3, bwR, bwG, bwB);
figure, imshowpair(img, overlay, 'montage');
title('Final Rendering of Centromere and Foci on SC')

% disp('Pausing for 4 seconds before closing images. Click on this window and hit ctrl-C to preserve them')
% pause(4)
% close all

%% Aberrant Detection - Minimum Eigenvalue Method
%Determine length of chromosome and distance from centromere to all foci

numCorners = zeros(1,redFound.NumObjects);
areas = zeros(1, redFound.NumObjects); % for next step
for i = 1:redFound.NumObjects % isolates each chromosome found
    blank = logical(a); % creates blank logical matrix
    blank(redFound.PixelIdxList{i}) = 1; % assigns all the listed pixels to 1
    
    % Crop and center
    cropped = regionprops(blank, 'image'); % tight crop
    cropped = cropped.Image;
    cropped = imtranslate(cropped,[20, 20],'OutputView','full'); % expands
    cropped = imtranslate(cropped, [-10,-10]); % centers
    
    % smooth it so smaller corners aren't detected?
    smooth = imgaussfilt(double(cropped), 1);
    
    % Detect Corners
    minEigen = 0.4275; % default...
    missing = 0;
    if redFound.NumObjects < 20 && redFound.NumObjects >= 12 % be more sensitive if missing any
        missing = 0.03*(20 - redFound.NumObjects);
    elseif redFound.NumObjects < 12 % can't let minEigen go too low
        minEigen = 0.1;
    end
    
    % Corner detection should become progressively more sensitive the more
    % missing chromosomes there are.
    corners = detectMinEigenFeatures(smooth, 'minQuality', minEigen - missing);
    
    numCorners(i) = corners.length(); % 5 seems like a good cut off
    areas(i) = size(redFound.PixelIdxList{i}, 1); % for next step
end

corner_deviants = find(numCorners >= 5); % empirically chose 5 @ minEigen = 0.4275
display(corner_deviants)

figure
corners = logical(a); % creates blank logical matrix
for x = 1:length(corner_deviants) %cycles thru aberrants
    corners(redFound.PixelIdxList{corner_deviants(x)}) = 1;
    imshowpair(img, corners, 'montage'),
    title('Potential Aberrants, Corners')
end

%% Aberrant Detection - Area Method
% create graph showing distribution
figure, bar(sort(areas));
refline(0,median(areas) - (1 - 0.075*missing)*iqr(areas)) %less sensitive
refline(0,median(areas))
refline(0,median(areas) + 0.75*iqr(areas))

cutoffLow = median(areas) - (1 - 0.075*missing)*iqr(areas); %less sensitive
cutoffHigh = median(areas) + 0.75*iqr(areas);
area_deviants = find(areas <= cutoffLow | areas >= cutoffHigh);
display(area_deviants)

figure
areaList = logical(a); % creates blank logical matrix
for x = 1:length(area_deviants) %cycles thru aberrants
    areaList(redFound.PixelIdxList{area_deviants(x)}) = 1;
    imshowpair(img, areaList, 'montage'),
    title('Potential Aberrants, Area')
end
%% Measurements?
figure
measures = zeros(1,redFound.NumObjects);
for i = 1:redFound.NumObjects % isolates each chromosome found
    blank = logical(a); % creates blank logical matrix
    blank(redFound.PixelIdxList{i}) = 1; % assigns all the listed pixels to 1
    
    % use perimeter function to get estimate
    measure = regionprops(blank, 'perimeter');
    perimFunction = measure.Perimeter/2;
    
    % outline
    [~, threshold] = edge(blank, 'sobel');
    perim = edge(blank,'sobel', threshold * fudgeFactor);
    
    imshow(perim)
    perimSum = sum(sum(perim))/2;
    %     display(perimFunction)
    %     display(perimSum)
    % Which method is more accurate? rip.
    measures(i) = (perimSum - perimFunction);
    
end
disp(measures)

%% TODO
%
% Detect intersections and overlap between chromosomes
%   - aberrant detection for number of centromeres per PixelxIDList blob
%   - use the area to decide if the threshold for blue channel needs to be bumped
%         - make the above work WITH the blue channel's new overlay
%         reduction
% MEASURING LENGTH
%   - Perimeter/2 - the initial curve of each end
%       - evaluate the percent error of the diagonal vs horizontal issue by
%       using snipped yarn on solid black background
%
%   x implement light centromere background reduction of (bwR & bwB)
%   x implement area approach
%   x Harris Corner Detector/Minimum Eigen Value
%         - run on each isolated SC, count total num
%   x Adjust image for the pixel intensity gradient (signal strength gets
%       weaker from left to right - adjust)?
%   ? attempt hough voting approach on the purely thresholded image
%       -low priority because corner approach better
%   ? imageJ macro call?
%       -not expecting high pay off
%
%
% PROBLEMS:
% - overlap
%     - "how many here" 3 -> 3 new additions to pixelID list
%     - using interaction with graphed data to select which pixelID to
%     delete (afterward)
% - touching
% - disconnected chromosome
% - unfilled in (partial as well)
%     - smaller area, longer length
%     - isolate and dilate ("erosion" techniques)
%
% - Using regionprops
%     - 'image'
%     - 'centroid'
%     - 'perimeter'
%
% -Commit to Git!
%
% https://www.mathworks.com/help/matlab/data_analysis/interacting-with-graphed-data.html
%
%
%   ? GUI
%          - XY, Aberrant, Accept?
%                 -Aberrant -> Overlap, touching, disconnected, unfilled in
%                     -when drawing overlap, allow to designate as XY
