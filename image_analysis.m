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

% pause(4)
% close

%% Detection of Red
% enhance the red of the image - counter the gradual decrease in red
% intensity from left to right
stretch = decorrstretch(red,'tol',0.02);
darker = imadjust(stretch, stretchlim(stretch),[0.9 0.99]);

bw = imbinarize(darker, graythresh(darker)); % binarizes with best threshold
bwR = bwareaopen(bw, 50); % reduces background noise of binarized image

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

% pause(2)
% close
%% Detection of blue
% females have (19 autosomal + 1 sex) centromeres == 20
% males have (19 autosomal + 2 WEAK) centromeres == 19

if ~flag
    numOfCentromeres = 20 - isMale; % searches for 20 if F, 19 if M
else % if bad image flag has been set, assume 1 chromosome = 1 centromere
    numOfCentromeres = redFound.NumObjects;
end

%defaults
numFound = 0;
threshold = graythresh(blue);
bckgrndReduct = 16; % to be adjusted

while numFound ~= numOfCentromeres
    bw = imbinarize(blue, threshold); % binarizes with best threshold
    bwB = bwareaopen(bw, bckgrndReduct); % reduces background noise of binarized image
    
    blueFound = bwconncomp(bwB, 8); % attempts to count number of objects
    numFound = blueFound.NumObjects;
    
    %adjust
    if numFound < numOfCentromeres && threshold > 0.005 % threshold is too high... not enough
        threshold = threshold - 0.005;
    elseif numFound > numOfCentromeres && threshold < 0.995% threshold is too low, too many
        threshold = threshold + 0.005;
    elseif threshold >= 0.995 %  max threshold used - adjust bckgrnd reduct.
        warning('Exceptional Behavior: Max threshold used - Increasing background reduction')
        bckgrndReduct = bckgrndReduct + 1;
        threshold = graythresh(blue); %reset
    end
    
    if bckgrndReduct > 30
        warning('Repeated Exceptional Behavior - Continuing with best image')
        break
    end
end

figure, imshowpair(img, label2rgb(labelmatrix(blueFound)), 'montage')
title(strcat(['Binarized Blue Channel, Centromeres identified = ', num2str(numFound),...
    ', threshold = ', num2str(threshold)]))

% pause(2)
% close

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
bw = imdilate(bw, [se90 se0]);

% prune away foci that don't overlap red
bwG = bwareaopen((bwR & bw), 8);

greenFound = bwconncomp(bwG, 8); % attempts to count number of objects
numFound = greenFound.NumObjects;

figure, imshowpair(img, bwG, 'montage'),
title(strcat(['Binarized Green Channel, Foci identified = ', num2str(numFound)]))
% pause(4)
% close all
%% Review the New Composite
% Simplify the image by only green and blue areas that overlap red
overlay = cat(3, bwR, (bwR & bwG), (bwB & bwR ));
figure, imshowpair(img, overlay, 'montage');
title('Final Rendering of Centromere and Foci on SC')

if flag
    error('A bad image was given. All windows will be left open. Ending.')
end

disp('Pausing for 4 seconds before closing images. Click on this window and hit ctrl-C to preserve them')
pause(4)
close all
%% Measurements
%Determine length of chromosome and distance from centromere to all foci

for i = 1:redFound.NumObjects % isolates each chromosome found
    
    blank = logical(a); % creates blank logical matrix
    blank(redFound.PixelIdxList{i}) = 1; % assigns all the listed pixels to 1
    skeleton = bwmorph(blank,'skel',Inf); % skeletonizes
    
    % get coordinates of all points on the skeleton - use to polyfit
    [x,y] = meshgrid(1:size(skeleton,1), 1:size(skeleton,2));
    result = [x(:),y(:),skeleton(:)]; % create list of x and y coordinates + binary value
    coords = result(result(:,3) == 1, [1 2]); % create list of all 1's
    
    % determine best polynomial representation
    %     R = 0; % https://www.mathworks.com/matlabcentral/newsreader/view_thread/114373
    aggr = zeros(1,19);
    for n = 2:20 % brute force test polynomial degrees 2 thru 20
        [p, stretch] = polyfit(coords(:,1),coords(:,2),n); % x, y , degree
        R = 1 - stretch.normr^2/norm(coords(:,2)-mean(coords(:,2)))^2; % the R squared value
        aggr(n-1) = R; % compile list of R-squared values
    end
    
    [v, degree] = max(aggr); % help? This is the best polynomial + its R-squared value
    fprintf('The Polynomial that best describes chromosome %i is degree %i and has an R-squared value of %0.5f\n', i, degree, R)
    
    imshow(skeleton);
    pause(0.5)
end



% TODO
%
% Detect intersections and overlap between chromosomes
%   - imageJ macro call?
%
% https://www.mathworks.com/help/matlab/data_analysis/interacting-with-graphed-data.html
% 
% Adjust image for the pixel intensity gradient (signal strength gets
% weaker from left to right - adjust)?
%