%% Image Analysis
% location = input('Supply (in single quotes) the filepath to image folder');
% img = imread(location);
close all
location = '~/Desktop/College/Research/PayseurLab/female.tif'; % DELETE
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
% montage(images,'Size', [1 1]), title('Raw Color Channels');


%% Detection of Red

% normalize strength of red channel (attempt to)
stretch = decorrstretch(red,'tol',0.02);
darker = imadjust(stretch, stretchlim(stretch),[0.02 0.99]);

% clean image
bw = imbinarize(darker, graythresh(darker)); % binarizes with best threshold
bwSharp = bwareaopen(bw, 50); % reduces background noise of binarized image
bwSmooth = imgaussfilt(double(bwSharp), 1); % makes the image smooth so that the fill method doesn't get confused
bwR = imbinarize(bwSmooth, graythresh(darker)); % binarizes again

% evaluate
redFound = bwconncomp(bwR, 8); % make a preliminary count of the objects
numFound = redFound.NumObjects;

% adjust image if needed
% https://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
[~, threshold] = edge(bwR, 'sobel');
dilationFactor = 90;
fudgeFactor = .5;
dilationOn = 3; %little difference btwn 2 and 3, whereas 1 == off

BWs = edge(bwR,'sobel', threshold * fudgeFactor); % detect edges
while numFound ~= 20
    seBegin = strel('line', dilationOn, dilationFactor); % may combine very close chromosomes
    seEnd = strel('line', dilationOn, 0); % but also will combine fragmented ones
    BWsdil = imdilate(BWs, [seBegin seEnd]); % dilates edges
    final = imfill(BWsdil, 'holes'); % fill in area within edges
    
    adjusted = bwconncomp(final, 8); % attempts to count number of objects
    numFound = adjusted.NumObjects;
    
    % attempt to adjust
    if numFound < 20 && dilationFactor > 0 % too much dilation
        dilationFactor = dilationFactor - 5;
        continue
    elseif numFound > 20 && dilationFactor < 120 % too little dilation
        dilationFactor = dilationFactor + 5;
        continue
    elseif numFound < 20 && dilationFactor == 0 && dilationOn ~= 1 % turns off dilation completly
        dilationOn = 1;
        continue
    elseif numFound > 20 && dilationFactor == 120 && dilationOn < 5
        dilationFactor = 50;
        dilationOn = dilationOn + 1;
        continue
    elseif (numFound < 20 && dilationOn == 1) || dilationFactor >= 120
        warning('Bad image')
        beep
        flag = true;
        break
    end
    
    % this statement will only be reached when loop is exiting
    % if dilation fixed anything, uses the adjusted image
    if(adjusted.NumObjects ~= redFound.NumObjects)
        redFound = adjusted;
        fprintf('dilation used')
    end
end

% figure, imshowpair(img, label2rgb(labelmatrix(redFound)), 'montage')
% title(strcat(['Binarized Red Channel, Chromosomes identified = ', num2str(numFound)]))

%% Detection of blue
% females have (19 autosomal + 1 sex) centromeres == 20
% males have (19 autosomal + 2 WEAK) centromeres == 19

if ~flag || redFound.NumObjects < 20 % second condition triggered by lots of overlaps
    numOfCentromeres = 20 - isMale; % searches for 20 if F, 19 if M
else % if bad image flag has been set, assume 1 chromosome = 1 centromere
    numOfCentromeres = redFound.NumObjects;
end

% defaults
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
    
    % Eliminate background noise and coincidental overlap of blue and red
    bwBRaw = bwareaopen(bw, bckgrndReduct);
    bwB = bwareaopen(bwBRaw & bwR, bckgrndReduct);
    
    % count number of centromeres found
    blueFound = bwconncomp(bwB, 8); % attempts to count number of objects
    numFound = blueFound.NumObjects;
    
    % adjust
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

% figure, imshowpair(img, label2rgb(labelmatrix(blueFound)), 'montage')
% title(strcat(['Binarized Blue Channel, Centromeres identified = ', ...
%     num2str(numFound), ', threshold = ', num2str(threshold)]))

%% Detection of green
% binarize image with best threshold
bw = imbinarize(green, graythresh(green));

% Special case: The foci are so faint that a bad threshold is selected
percent_SC_covered = sum(sum(bw & bwR))/ sum(sum(bwR)); % if too many foci
if(percent_SC_covered > 0.3)
    warning('special case thresholding: green')
    bw = imbinarize(green, 0.175); % arbitrary threshold
end

%reduce background noise
bw = bwareaopen(bw, 4);

% find foci
[~, threshold] = edge(bw, 'sobel');
bw = edge(bw,'sobel', threshold * fudgeFactor);

% enlarge
seBegin = strel('line', 2, 90);
seEnd = strel('line', 2, 0);
lines = imdilate(bw, [seBegin seEnd]);
bw = imfill(lines, 'holes');

% prune away foci that don't overlap red
bwG = bwareaopen((bwR & bw), 8);

% attempt to count number of foci
greenFound = bwconncomp(bwG, 8);
numFound = greenFound.NumObjects;

% figure, imshowpair(img, bwG, 'montage'),
% title(strcat(['Binarized Green Channel, Foci identified = ', num2str(numFound)]))

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

corners_list = logical(a); % creates blank logical matrix
for i = 1:length(corner_deviants) %cycles thru aberrants
    corners_list(redFound.PixelIdxList{corner_deviants(i)}) = 1;
end

% figure, imshow(corners_list),title('Potential Aberrants, Corners')
%% Aberrant Detection - Area Method
% create graph showing distribution
% figure, bar(sort(areas));
% refline(0,median(areas) - (1 - 0.075*missing)*iqr(areas)) %less sensitive
% refline(0,median(areas))
% refline(0,median(areas) + 0.75*iqr(areas))

cutoffLow = median(areas) - (1 - 0.075*missing)*iqr(areas); %less sensitive
cutoffHigh = median(areas) + 0.75*iqr(areas);
area_deviants = find(areas <= cutoffLow | areas >= cutoffHigh);

area_list = logical(a); % creates blank logical matrix
for i = 1:length(area_deviants) %cycles thru aberrants
    area_list(redFound.PixelIdxList{area_deviants(i)}) = 1;
end

% figure, imshow(area_list),title('Potential Aberrants, Area')
%% Aberrant Detection - Centromere method
centromere_deviants = zeros(1,redFound.NumObjects);
for i = 1:redFound.NumObjects
    blank = logical(a); % creates blank logical matrix
    blank(redFound.PixelIdxList{i}) = 1; % plots chromosome
    centromeres = bwconncomp(bwB & blank, 8); % plots only the overlap
    
    %evaluate if aberrant
    if(centromeres.NumObjects ~= 1) 
        centromere_deviants(i) = i;
    end
end
centromere_deviants = find(centromere_deviants ~= 0); % clean out zeros

centromere_list = logical(a);
for i = 1:length(centromere_deviants) 
    centromere_list(redFound.PixelIdxList{centromere_deviants(i)}) = 1;
end
% figure,imshow(centromere_list),title('Potential Aberrants, Centromere')

%% User Input?
all_deviants = unique(horzcat(area_deviants, corner_deviants,centromere_deviants)); % collects aberrants

aberrants = logical(a); % creates blank logical matrix
for i = 1:length(all_deviants) % cycles thru known aberrants
    aberrants(redFound.PixelIdxList{all_deviants(i)}) = 1; % marks them on plot
end

aberrants = imcomplement(aberrants); % flips colors so easier to use cross-hair
refresh = aberrants;

clear('true_aberrants') % clears out any previous data
done = false; % sentinel value for completion
while(~done)
    
    aberrants = refresh; % makes sure the user receives a fresh image
    
    % customizes display
    beep
    figure, imshow(aberrants), title('Please click each aberrant silhouette. Use the delete key if you make any mistakes - it will remove the most recent click. Hit enter when done.')
    figh = gcf;
    pos = get(figh,'position');
    set(figh,'position',[pos(1:2)*0.5 pos(3:4)*3]);
    
    % receive user input
    [x,y, buttons] = ginput;
    user_input = uint64([x,y]); % cast this because it needs to be an integer
    
    % if the user hits the delete key
    if(sum(ismember(buttons, 8)))
        toRemove = find(buttons == 8) - 1; % remove the click before
        
        if(~sum(ismember(toRemove, 0))) % checks for delete key as first press
            user_input(toRemove,:) = 0; % marks the bad click
            user_input(toRemove + 1,:) = 0; % marks the delete-key press
            user_input(ismember(user_input, [0,0],'rows'),:) = []; % deletes
            % corrects coordinates too
            x = user_input(:,1);
            y = user_input(:,2);
        else
            uiwait(msgbox('Delete-key pressed too early.','Warning','modal'));
            close(gcf)
            continue
        end
    end
    
    % Create circles where the user clicks
    if(size(user_input,1) > 1)
        for i = 1:length(user_input) % marks user input
            aberrants = insertShape(double(aberrants),'circle',[x(i),y(i),3], 'color', 'red');
        end
    elseif size(user_input,1) == 1
        aberrants = insertShape(double(aberrants),'circle',[x,y,3], 'color', 'red');
    elseif size(user_input,1) == 0
        break % user didn't want to mark any aberrants
    end
    
    % See which options the user clicked on
    true_aberrants = zeros(1,redFound.NumObjects); % pre-allocate array for matches
    for i = 1:length(all_deviants) % see what objects the user clicked on
        
        [y,x] = ind2sub(size(a),redFound.PixelIdxList{all_deviants(i)}); % somehow this returns y then x....
        coords = [x,y];
        
        match = intersect(user_input,coords,'rows'); % do any coordinates on this silhouette match user input?
        if ~isempty(match)
            true_aberrants(i) = all_deviants(i); % adds confirmed aberrant to list
            aberrants = insertShape(double(aberrants),'circle',[match(1,1),match(1,2),3],...
                'color', 'green','LineWidth', 2); % confirms it on picture
            
            if(size(match,1) > 1) % if user clicked same chromosome multiple times, turns those green
                for j = 2:(size(match,1))
                    aberrants = insertShape(double(aberrants),'circle',[match(j,1),match(j,2),3],...
                        'color', 'green','LineWidth', 2); % confirms it on picture
                end
            end
            
            % Only keep track of user input that doesn't find a match
            toRemove = find(ismember(user_input, match,'rows'));
            if(~isempty(toRemove)) % in case user clicks on same chromosome twice
                user_input(toRemove, :) = [];
            end
        end
    end
    true_aberrants(true_aberrants==0) = []; % remove all 0's
    hold on, imshow(aberrants), hold off % show the user's clicks
    
    % invalid clicks are not allowed, restart if there are any
    if(~isempty(user_input))
        uiwait(msgbox('Not all user input points were found.','Warning','modal'));
        close(gcf)
    else
        done = true;
    end
end


hold on, imshow(aberrants), hold off % show on same plot
if(exist('true_aberrants','var'))
    fprintf('you have selected chromosome %i\n', true_aberrants);
else
    true_aberrants = 0; % initialize
end
pause(1)
close(gcf)
%% Measurements?
outlines = cell(1, redFound.NumObjects);
image_data = struct('Chromosome_Length', {zeros(redFound.NumObjects,1)},...
    'Foci_Distances', {[]});
for i = 1:redFound.NumObjects % isolates each chromosome found
    blank = logical(a); % creates blank logical matrix
    blank(redFound.PixelIdxList{i}) = 1; % assigns all the listed pixels to 1
    
    % use perimeter function to get estimate
    measure = regionprops(blank, 'perimeter');
    perimFunction = measure.Perimeter/2;
    
    % outline
    [~, threshold] = edge(blank, 'sobel');
    perim = edge(blank,'sobel', threshold * fudgeFactor);
    outlines(i) = {double(perim)}; % save outline for later use
    
    perimSum = sum(sum(perim))/2;
    image_data.Chromosome_Length(i) = perimSum;
end
%% Dijkstra's Approach to automated measuring of inter-foci distance

% create a map to pass Dijkstra's method
overall = overlay; % copies the type of variable
overall(:) = 0; % clears it out
for i = 1:redFound.NumObjects
    if(sum(true_aberrants == i)) % will skip any aberrants
        image_data.Foci_Distances{i} = [];
        continue
    end
    
    body = logical(a); % blank map of 0's
    body(redFound.PixelIdxList{i}) = 1; % body of chromosome
    map = 2*double(body); % begins map creation
    skeleton = double(bwmorph(map,'skel',Inf)); % saves skeletal chromosome
    
    map = map + 2*outlines{i}; % makes the perimeter higher cost
    map(skeleton == 1) = 1; % makes the skeletal path the lowest cost
    
    % display - DELETE
    blank = logical(a); blank(redFound.PixelIdxList{i}) = 1;
    overall = overall + cat(3, blank, logical(skeleton), logical(outlines{i}));
    imshowpair(img, overall, 'montage'); hold on
    
    % find centroid of centromere
    centromere = regionprops(body & bwB, 'centroid');
    if(isempty(centromere)) % skip this chromosome if there aren't any centromeres
        warning('ABERRANT GOT THROUGH')
        continue,
    end
    centromere = struct2cell(centromere);
    centromere = int64(centromere{1}); % unpack the structure
    %     overall = insertShape(overall,'circle',[centromere(1),centromere(2),6], 'color', 'blue'); % DELETE
    
    % find centroids of foci
    foci = regionprops(body & bwG, 'centroid');
    if(isempty(foci)), continue, end % skip this chromosome if there aren't any foci
    foci = int64(cat(1, foci.Centroid));
    %     overall = insertShape(overall,'circle',[foci(:,1),foci(:,2),repmat(6,size(foci,1),1)], 'color', 'green'); % DELETE
    %     imshow(overall)
    
    [distance, path] = dijkstra_image(map, centromere, foci);
    image_data.Foci_Distances{i} = distance;
end
if(size(image_data.Foci_Distances, 1) == 1)
    image_data.Foci_Distances = image_data.Foci_Distances'; % makes data organization consistent
end

return
% body = double(body);
% body(path) = 0;
% body = imcomplement(body);
% body = insertShape(body,'circle',[foci(:,1),foci(:,2),repmat(3,size(foci,1),1)], 'color', 'green'); % DELETE
% body = insertShape(body,'circle',[centromere(1),centromere(2),3], 'color', 'blue'); % DELETE
% imshow(body) % DELETE
%% Spline Measuring
for i = 1:length(true_aberrants)
    blank = logical(a); % creates blank logical matrix
    blank(redFound.PixelIdxList{true_aberrants(i)}) = 1; % assigns all the listed pixels to 1
    
    % Crop and center
    cropped = regionprops(blank, 'image'); % tight crop
    cropped = cropped.Image;
    cropped = imtranslate(cropped,[20, 20],'OutputView','full'); % expands
    cropped = imtranslate(cropped, [-10,-10]); % centers
    cropped = imcomplement(cropped);
    imshow(cropped)
    
    figh = gcf;
    pos = get(figh,'position');
    set(figh,'position',[pos(1:2)*0.5 pos(3:4)*3]);
    
    % splining
    [x,y] = ginput;
    X = [x,y];
    for j = 2:size(X,1)
        if j == 2, distance = 0; end
        distance = distance + pdist([X(j-1,:);X(j,:)], 'euclidean');
    end
    display(distance)
end
% TODO: store distances
%using splines to measure aberrants
%evaluate two objects of same length, diagonal and horizontal

%% TODO
%
%   - Detect intersections and overlap between chromosomes
%       - aberrant detection for number of centromeres per PixelxIDList blob
%       - use the area to decide if the threshold for blue channel needs to be bumped
%       - make the above work WITH the blue channel's new overlay
%         reduction
%
%   - evaluate the percent error of the diagonal vs horizontal issue by
%       using snipped yarn on solid black background
%
% Problem Cases:
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
%   ? GUI
%          - XY, Aberrant, Accept, Delete?
%                 -Aberrant -> Overlap, touching, disconnected, unfilled in
%                     -when drawing overlap, allow to designate as XY
