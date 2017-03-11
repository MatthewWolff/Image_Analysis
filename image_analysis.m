%% Image Analysis
% location = input('Supply (in single quotes) filepath to an image');
% img = imread(location);
location = '~/Desktop/College/Research/PayseurLab/male.tif'; % DELETE
flag = false; %DELETE - to mark when a bad image has been given
% in the case of a bad image, analysis will continue in spite of that by
% disregarding the previous expected values (e.g. centromere #)

img = imread(location);

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
dilationOn = 3; %little differnce btwn 2 and 3, whereas 1 == off

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
bckgrndReduct = 16; % to be adjusted
adjustCount = 0;
adjustLimit = 2000; % to catch infinite loops -> normal is < 100

while numFound ~= numOfCentromeres
    bw = imbinarize(blue, threshold); % binarizes with best threshold
    bwB = bwareaopen(bw, bckgrndReduct); % reduces background noise of binarized image
    
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
bw = bwareaopen(bw, 4); % reduces background noise of binarized image

% find foci
[~, threshold] = edge(bw, 'sobel');
bw = edge(bw,'sobel', threshold * fudgeFactor);

% enlarge
se90 = strel('line', 2, 45);
se0 = strel('line', 2, 0);
lines = imdilate(bw, [se90 se0]);
bw = imfill(lines, 'holes');

% prune away foci that don't overlap red
bwG = bwareaopen((bwR & bw), 8);

greenFound = bwconncomp(bwG, 8); % attempts to count number of objects
numFound = greenFound.NumObjects;

figure, imshowpair(img, bwG, 'montage'),
title(strcat(['Binarized Green Channel, Foci identified = ', num2str(numFound)]))

%% Review the New Composite
% Simplify the image by only green and blue areas that overlap red
overlay = cat(3, bwR, (bwR & bwG), (bwB & bwR ));
figure, imshowpair(img, overlay, 'montage');
title('Final Rendering of Centromere and Foci on SC')

% if flag
%     error('A bad image was given. All windows will be left open. Ending.')
% end

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
    skeleton = bwmorph(blank,'skel',Inf); % skeletonizes
    
    % get coordinates of all points on the skeleton DELETE - unused
    [x,y] = meshgrid(1:size(skeleton,1), 1:size(skeleton,2));
    result = [x(:),y(:),skeleton(:)]; % create list of x and y coordinates + binary value
    coords = result(result(:,3) == 1, [1 2]); % create list of all 1's
    
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
        missing = 0.05*(20 - redFound.NumObjects);
    elseif redFound.NumObjects < 12 % can't let minEigen go too low
        minEigen = 0.1;
    end
    
    % Corner detection should become progressively more sensitive the more
    % missing chromosomes there are.
    corners = detectMinEigenFeatures(smooth, 'minQuality', minEigen - missing);
    
    numCorners(i) = corners.length(); % 5 seems like a good cut off
    areas(i) = size(redFound.PixelIdxList{i}, 1); % for next step
    
    %     imshow(smooth); hold on;
    %     plot(corners.selectStrongest(4));
    %     hold off;
    
    
end

% corners = max(numCorners) % only finds the first max
corner_deviants = find(numCorners >= 5); % empirically chose 5 @ minEigen = 0.4275
display(corner_deviants)

corners = logical(a); % creates blank logical matrix
for x = 1:length(corner_deviants) %cycles thru aberrants
    corners(redFound.PixelIdxList{corner_deviants(x)}) = 1;
    imshowpair(img, corners, 'montage'),
    title('Potential Aberrants')
end

%% Aberrant Detection - Area Method
close all

% create graph showing distribution
figure, bar(sort(areas));
refline(0,median(areas) - iqr(areas))
refline(0,median(areas))
refline(0,median(areas) + iqr(areas))

cutoffLow = median(areas) - iqr(areas);
cutoffHigh = median(areas) + iqr(areas);
area_deviants = find(areas <= cutoffLow | areas >= cutoffHigh);

figure
areaList = logical(a); % creates blank logical matrix
for x = 1:length(area_deviants) %cycles thru aberrants
    areaList(redFound.PixelIdxList{area_deviants(x)}) = 1;
    imshowpair(img, areaList, 'montage'),
    title('Potential Aberrants')
end


% TODO
%
% Detect intersections and overlap between chromosomes
%   - implement area approach
%   - check number of centromeres per pixelID obj
%   x Harris Corner Detector
%         - run on each isolated SC, count total num
%   x Adjust image for the pixel intensity gradient (signal strength gets
%       weaker from left to right - adjust)?
%   ? attempt hough voting approach on the purely thresholded image
%       -low priority because corner approach better
%   ? imageJ macro call?
%       -not expecting high pay off
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
%