function [extracted_data,image_data] = image_analysis(location)

% % Image Analysis
clearvars -except location
close all
% location = '~/Desktop/College/Research/PayseurLab/female.tif'; % DELETE
% location = 'C:\Users\alpeterson7\Documents\matt wolf\WSB1 (5).tif'; % DELETE

img = imread(location);
flag = false; % DELETE - to mark when a bad image has been given
% in the case of a bad image, analysis will continue in spite of that by
% disregarding the previous expected values (e.g. centromere #)

% Determine if the slide is male or female from filepath
[~,name,~] = fileparts(location);
isFemale = regexp(name,'female','ONCE');
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

% images = [img, just_red; just_green, just_blue];
% montage(images,'Size', [1 1]), title('Raw Color Channels');
resize_window = @(positions, magnitude) [positions(1:2)-50*magnitude, positions(3:4)+100*magnitude];

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
        %         warning('Bad image')
        %         beep
        flag = true;
        break
    end
    
    % this statement will only be reached when loop is exiting
    % if dilation fixed anything, uses the adjusted image
    if(adjusted.NumObjects ~= redFound.NumObjects)
        redFound = adjusted;
        warning('dilation used')
    end
end
%
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
blueFlag = false;
while numFound ~= numOfCentromeres
    % Prepare Image for counting
    bw = imbinarize(blue, threshold); % binarizes with best threshold
    
    %Special Case: Auto-thresholding was very bad - arbitrarily make 0.7
    percentOfImg = sum(sum(bw))/ (size(img, 1) * size(img, 2));
    if(adjustCount == 0 && percentOfImg > 0.2)
        warning('special case thresholding: blue')
        threshold = 0.7;
        blueFlag = true;
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
        if(blueFlag)
            threshold = 0.7;
        else
            threshold = graythresh(blue); %reset
        end
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

% Special case: There's so much noise around the cell that we're finding
% too many foci
numFound = 100;
greenthresh = 0.25;
while(numFound > 40 && greenthresh < 1)
    
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
    
    if(numFound > 40)
        bw = imbinarize(green, greenthresh);
        greenthresh = greenthresh + 0.05;
    end
    
end
% figure, imshowpair(img, bwG, 'montage'),
% title(strcat(['Binarized Green Channel, Foci identified = ', num2str(numFound)]))

%% Review the New Composite
% Simplify the image by only green and blue areas that overlap red
overlay = double(cat(3, bwR, bwG, bwB));
figure, imshowpair(overlay, img, 'montage');
title('Final Rendering of Centromere and Foci on SC')

figh = gcf;
pos = get(figh,'position');
set(figh,'position',resize_window(pos,6));

% Construct a questdlg with three options
choice = questdlg('Would you like to add foci?', ...
    'Optional Foci Addition','Yes','No','No'); % last option is default

refresh = imcomplement(overlay); % for visual ease

done = false; % sentinel value for completion
if(strcmp(choice,'No')) % skip if user doesn't want to add foci
    done = true;
end
while(~done)
    
    temp = refresh; % makes sure the user receives a fresh image
    
    % customizes display
    beep
    figure, imshow(temp), title('Please click where you would like to add foci. Press enter when done.')
    figh = gcf;
    pos = get(figh,'position');
    set(figh,'position',resize_window(pos,6));
    
    % receive user input
    [x,y, buttons] = ginput;
    user_input = uint32([x,y]); % cast this because it needs to be an integer
    
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
            temp = insertShape(double(temp),'circle',[x(i),y(i),3], 'color', 'black');
        end
    elseif size(user_input,1) == 1
        temp = insertShape(double(temp),'circle',[x,y,3], 'color', 'black');
    elseif size(user_input,1) == 0
        break % user didn't want to mark any aberrants
    end
    
    hold on, imshow(temp), hold off % show the user's clicks
    
    
    choice = questdlg('Add foci to these locations?', ...
        'Foci Confirmation','Yes','Redo','Yes'); % last option is default
    
    % invalid clicks are not allowed, restart if there are any
    if(strcmp(choice,'Yes'))
        done = true;
    else
        close(gcf)
        continue; %re-do
    end
    
    % create the foci
    new_foci = logical(a);
    x = user_input(:,1); y = user_input(:,2);
    point = [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y+1;x+1,y+1;x+1,y;x+1,y-1];
    center = sub2ind(size(new_foci), point(:,2),point(:,1));
    new_foci(center) = 1;
    
    for iterations = 1:4 % grow the centromere
        % prepare to dilate
        [~, threshold] = edge(new_foci, 'sobel');
        new_foci = edge(new_foci,'sobel', threshold * fudgeFactor);
        
        % dilate and fill the drawn line
        seBegin = strel('line', 2, 100);
        seEnd = strel('line', 2, -1);
        center = imdilate(new_foci, [seBegin seEnd]);
        new_foci = imfill(center, 'holes');
    end
    
    bwG(new_foci) = 1; % add to pre-existing foci
    bwG = bwG & bwR; % clean up (remove any green that doesn't overlap red)
    
    overlay = double(cat(3, bwR, bwG, bwB)); % update
end
close(gcf)

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
    missing = (20 - redFound.NumObjects);
    missing_adjustment = 0;
    if redFound.NumObjects < 20 && redFound.NumObjects >= 12 % be more sensitive if missing any
        missing_adjustment = abs(0.03*missing);
    elseif redFound.NumObjects < 12 % can't let minEigen go too low
        minEigen = 0.1;
    end
    
    % Corner detection should become progressively more sensitive the more
    % missing chromosomes there are.
    corners = detectMinEigenFeatures(smooth, 'minQuality', minEigen - missing_adjustment);
    
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

%% Removal of Abnormally Large Centromeres
% detect abnormally large centromeres
centromere_deviants = zeros(1,redFound.NumObjects); % initialize
centromere_sizes = cellfun('length',blueFound.PixelIdxList);
too_large = mean(centromere_sizes) + 2*iqr(centromere_sizes);
too_small = mean(centromere_sizes) - 1.5*iqr(centromere_sizes);
abnormal = find(centromere_sizes >= too_large | centromere_sizes <= too_small); % find offenders
if(~isempty(abnormal)) % if offenders present
    removed = 0;
    for i = 1:length(abnormal) % decide whether or not to remove them
        blank = logical(a);
        blank(blueFound.PixelIdxList{abnormal(i)}) = 1;
        combo = cat(3,imcomplement(bwR & blank),imcomplement(bwR),imcomplement(a));
        figure, imshowpair(combo,img,'montage')
        figh = gcf;
        pos = get(figh,'position');
        set(figh,'position',resize_window(pos,6));
        
        choice = questdlg('Remove abnormal centromere in blue?','Abnormal Centromere detected',...
            'Yes','Ignore','Ignore');
        if(strcmp(choice,'Yes'))
            % remove offender
            blueFound.PixelIdxList{abnormal(i)} = [];
            removed = removed + 1;
            
            % update blue channel to reflect removal
            new_blue = logical(a);
            for j = 1:blueFound.NumObjects
                new_blue(blueFound.PixelIdxList{j}) = 1;
            end
            bwB = new_blue;
        end
        close(gcf)
    end
    blueFound.NumObjects = blueFound.NumObjects - removed;
    blueFound.PixelIdxList = blueFound.PixelIdxList(~cellfun('isempty',blueFound.PixelIdxList));
end

%% Aberrant Detection - Centromere method
message_counter = 1; % set this to 0 after a message is given so it isn't repeated
for i = 1:redFound.NumObjects
    blank = logical(a); % creates blank logical matrix
    blank(redFound.PixelIdxList{i}) = 1; % plots chromosome
    centromeres = bwconncomp(bwB & blank, 8); % plots only the overlap
    
    %evaluate if aberrant
    if(centromeres.NumObjects == 1)
        continue;
    elseif(centromeres.NumObjects > 1)
        centromere_deviants(i) = i;
    end
    
    combo = cat(3,imcomplement(blank & bwB),imcomplement(blank),imcomplement(a));
    refresh = combo;
    
    beep
    blue_red = cat(3,red,a,blue);
    figure,imshowpair(combo, blue_red,'montage')
    if centromeres.NumObjects == 0
        title('A chromosome with 0 centromeres has been detected! Click on the image where the centromere should go, or hit enter to ignore.')
    else
        title('A chromosome with multiple centromeres has been detected! Hit Enter to ignore, or add one by clicking.')
    end
    
    myFig = gcf;
    pos = get(myFig,'position');
    set(myFig,'position',resize_window(pos,6));
    
    done = false; % sentinel value for completion
    while(~done)
        
        combo = refresh; % makes sure the user receives a fresh image
        hold on,imshowpair(combo, blue_red,'montage')
        if centromeres.NumObjects == 0
            title('A chromosome with 0 centromeres has been detected! Click on the image where the centromere should go, or hit enter to ignore.')
        else
            title('A chromosome with multiple centromeres has been detected! Hit Enter to ignore, or add one by clicking.')
        end
        myFig = gcf;
        set(myFig,'position',resize_window(pos,6));
        % customizes display
        if message_counter
            uiwait(msgbox('Please add/remove centromeres if & when necessary.','Instructions','modal'));
            message_counter = 0;
        end
        
        % receive user input
        [x,y, buttons] = ginput;
        user_input = uint32([x,y]); % cast this because it needs to be an integer
        
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
            for j = 1:length(user_input) % marks user input
                combo = insertShape(double(combo),'circle',[x(j),y(j),4], 'color', 'black');
            end
        elseif size(user_input,1) == 1
            combo = insertShape(double(combo),'circle',[x,y,4], 'color', 'black');
        elseif size(user_input,1) == 0
            centromere_deviants(i) = i;
            break % user didn't want to mark any aberrants
        end
        
        hold on, imshow(combo), hold off % show the user's clicks
        
        choice = questdlg('Add centromeres to these locations?', ...
            'Centromere Confirmation','Yes','Redo','Yes'); % last option is default
        
        if(strcmp(choice,'Yes'))
            done = true;
        else
            close(gcf)
            continue; %re-do
        end
        
        % create the centromere
        new_centromere = logical(a);
        x = user_input(:,2); y = user_input(:,1);
        point = [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y+1;x+1,y+1;x+1,y;x+1,y-1];
        center = sub2ind(size(new_centromere), point(:,1),point(:,2));
        new_centromere(center) = 1;
        
        for iterations = 1:7 % grow the centromere
            % prepare to dilate
            [~, threshold] = edge(new_centromere, 'sobel');
            new_centromere = edge(new_centromere,'sobel', threshold * fudgeFactor);
            
            % dilate and fill the drawn line
            seBegin = strel('line', 2, 100);
            seEnd = strel('line', 2, -1);
            center = imdilate(new_centromere, [seBegin seEnd]);
            new_centromere = imfill(center, 'holes');
        end
        bwB(new_centromere) = 1;
        bwB = bwB & bwR; % clean up
        image_data.Composite = cat(3,bwR,bwG,bwB);
    end
    
    hold on, imshowpair(combo, blue_red,'montage'), hold off
    pause(1)
    close(gcf)
    
end
centromere_deviants = find(centromere_deviants ~= 0); % clean out zeros

centromere_list = logical(a);
for i = 1:length(centromere_deviants)
    centromere_list(redFound.PixelIdxList{centromere_deviants(i)}) = 1;
end

% figure,imshow(centromere_list),title('Potential Aberrants, Centromere')

%% User Input - Aberrant Confirmation
all_deviants = unique(horzcat(area_deviants, corner_deviants,centromere_deviants)); % collects aberrants

aberrants = logical(a); % creates blank logical matrix
for i = 1:length(all_deviants) % cycles thru known aberrants
    aberrants(redFound.PixelIdxList{all_deviants(i)}) = 1; % marks them on plot
end

combo = cat(3,imcomplement(aberrants & bwB),imcomplement(aberrants),imcomplement(a));
refresh = combo;

clear('true_aberrants') % clears out any previous data
done = false; % sentinel value for completion
while(~done)
    
    aberrants = refresh; % makes sure the user receives a fresh image
    
    % customizes display
    beep
    figure, imshow(aberrants), title('Please click each aberrant silhouette. Use the delete key if you make any mistakes - it will remove the most recent click. Hit enter when done.')
    figh = gcf;
    pos = get(figh,'position');
    set(figh,'position',resize_window(pos,6));
    uiwait(msgbox('Please identify any aberrants.','Instructions','modal'));
    
    % receive user input
    [x,y, buttons] = ginput;
    user_input = uint32([x,y]); % cast this because it needs to be an integer
    
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
            
            % Only keeps track of user input that doesn't find a match
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
if(~exist('true_aberrants','var'))
    true_aberrants = []; % initialize
end
true_aberrants = true_aberrants'; % makes into a vertical vector
pause(1)
close(gcf)
%% Aberrant Classification
% create true list
close all
chromosomes = struct('PixelIdxList', {redFound.PixelIdxList},'NumObjects', redFound.NumObjects);
found_fragments = []; % keep track of which chromosomes are fragments
redraw = []; % keep track of which need to be redrawn
if(~isempty(true_aberrants)) % skip if there aren't aberrants
    figure % creates new window to use
    for i = 1:length(true_aberrants)
        blank = logical(a); % creates blank logical matrix
        blank(redFound.PixelIdxList{true_aberrants(i)}) = 1; % assigns all the listed pixels to 1
        figure(1),
        subplot(1,2,1),imshowpair(blank,bwR), title('Aberrant Designation')
        subplot(1,2,2),imshow(img), title('Original Image')
        
        figh = gcf;
        if(i == 1)
            pos = get(figh,'position');
            set(figh,'position',resize_window(pos,6));
            prompt = {'This box can be left open.';'';...
                'If the aberrant should be deleted, type "delete" or 0';...
                'If the aberrant simply needs to be redrawn, type 1';...
                'If there are multiple chromosomes, type the number (2 to 6)';...
                'If the aberrant is part of a new fragmented chromosome, type "b"';...
                'If the aberrant is part of known fragmented chromosome, type "c"'};
            
            if(isMale)
                prompt{length(prompt) + 1} = 'If the aberrant is the XY chromosome, type "xy"';
            else
                clear('XY') % if XY doesn't exist, Dijkstra's special case won't execute
            end
            msgbox(prompt)
        end
        
        valid = false;
        while(~valid)
            designation = inputdlg('Please select a designation for the aberrant (marked in white)');
            switch designation{1}
                case 'delete' % delete key - deletes entry
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects - 1;
                    valid = true;
                case '0' % 0.. same as delete
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects - 1;
                    valid = true;
                case '1' % 1
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    redraw = [redraw,true_aberrants(i)];
                    valid = true;
                case '2' % 2
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects + 1;
                    valid = true;
                    for r = 1:2, redraw = [redraw,true_aberrants(i)]; end
                case '3' % 3
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects + 2;
                    valid = true;
                    for r = 1:3, redraw = [redraw,true_aberrants(i)]; end
                case '4' % 4
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects + 3;
                    valid = true;
                    for r = 1:4, redraw = [redraw,true_aberrants(i)]; end
                case '5' % 5
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects + 4;
                    valid = true;
                    for r = 1:5, redraw = [redraw,true_aberrants(i)]; end
                case '6' % 6
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects + 5;
                    valid = true;
                    for r = 1:6, redraw = [redraw,true_aberrants(i)]; end
                case 'xy' % x - xy chromosome
                    XY = struct('Centromeres',zeros(2), 'PixelIdxList',...
                        redFound.PixelIdxList{true_aberrants(i)},'Length',-1);
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    
                    % temporarily shrink the list
                    chromosomes.NumObjects = chromosomes.NumObjects - 1;
                    
                    aberrants = logical(a); % creates blank logical matrix
                    aberrants(redFound.PixelIdxList{true_aberrants(i)}) = 1; % marks them on plot
                    
                    aberrants = imcomplement(aberrants); % flips colors so easier to use cross-hair
                    refresh = aberrants;
                    
                    done = false; % sentinel value for completion
                    while(~done)
                        
                        aberrants = refresh; % makes sure the user receives a fresh image
                        
                        % customizes display
                        beep
                        figure, imshow(aberrants), title('Please click on each end of the XY chromosome')
                        myFig = gcf;
                        pos = get(myFig,'position');
                        set(myFig,'position',resize_window(pos,6));
                        
                        % receive user input
                        [x,y, buttons] = ginput;
                        user_input = uint32([x,y]); % cast this because it needs to be an integer
                        
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
                            for j = 1:length(user_input) % marks user input
                                aberrants = insertShape(double(aberrants),'circle',[x((j)),y(j),3], 'color', 'red');
                            end
                        elseif size(user_input,1) == 1
                            aberrants = insertShape(double(aberrants),'circle',[x,y,3], 'color', 'red');
                        elseif size(user_input,1) == 0
                            break % user didn't want to mark any aberrants
                        end
                        
                        
                        
                        [y,x] = ind2sub(size(a),redFound.PixelIdxList{true_aberrants(i)}); % somehow this returns y then x....
                        coords = [x,y];
                        
                        match = intersect(user_input,coords,'rows'); % do any coordinates on this silhouette match user input?
                        if ~isempty(match)
                            XY.Centromeres = match;
                            aberrants = insertShape(double(aberrants),'circle',[match(1,1),match(1,2),3],...
                                'color', 'green','LineWidth', 2); % confirms it on picture
                            
                            if(size(match,1) > 1) % if user clicked same chromosome multiple times, turns those green
                                for k = 2:(size(match,1))
                                    aberrants = insertShape(double(aberrants),'circle',[match(k,1),match(k,2),3],...
                                        'color', 'green','LineWidth', 2); % confirms it on picture
                                end
                            end
                            
                            % Only keeps track of user input that doesn't find a match
                            toRemove = find(ismember(user_input, match,'rows'));
                            if(~isempty(toRemove)) % in case user clicks on same chromosome twice
                                user_input(toRemove, :) = [];
                            end
                        end
                        
                        hold on, imshow(aberrants), hold off % show the user's clicks
                        
                        % invalid clicks are not allowed, restart if there are any
                        if(~isempty(user_input))
                            uiwait(msgbox('Not all user input points were found.','Warning','modal'));
                            close(gcf)
                        else
                            done = true;
                        end
                    end
                    
                    hold on, imshow(aberrants), hold off
                    pause(1)
                    close(gcf)
                    valid = true;
                    
                case 'b' % b - broken chromosome
                    aberrants = logical(a); % creates blank logical matrix
                    for j = 1:length(all_deviants) % cycles thru known aberrants
                        aberrants(redFound.PixelIdxList{all_deviants(j)}) = 1; % marks them on plot
                    end
                    
                    aberrants = imcomplement(aberrants); % flips colors so easier to use cross-hair
                    refresh = aberrants;
                    
                    done = false; % sentinel value for completion
                    while(~done)
                        
                        aberrants = refresh; % makes sure the user receives a fresh image
                        
                        % customizes display
                        beep
                        figure, imshow(aberrants), title('Please click all fragments of the chromosome in question')
                        myFig = gcf;
                        pos = get(myFig,'position');
                        set(myFig,'position',resize_window(pos,6));
                        
                        % receive user input
                        [x,y, buttons] = ginput;
                        user_input = uint32([x,y]); % cast this because it needs to be an integer
                        
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
                            for j = 1:length(user_input) % marks user input
                                aberrants = insertShape(double(aberrants),'circle',[x((j)),y(j),3], 'color', 'red');
                            end
                        elseif size(user_input,1) == 1
                            aberrants = insertShape(double(aberrants),'circle',[x,y,3], 'color', 'red');
                        elseif size(user_input,1) == 0
                            break % user didn't want to mark any aberrants
                        end
                        
                        % See which options the user clicked on
                        fragments = zeros(1,redFound.NumObjects); % pre-allocate array for matches
                        for j = 1:length(all_deviants) % see what objects the user clicked on
                            
                            [y,x] = ind2sub(size(a),redFound.PixelIdxList{all_deviants(j)}); % somehow this returns y then x....
                            coords = [x,y];
                            
                            match = intersect(user_input,coords,'rows'); % do any coordinates on this silhouette match user input?
                            if ~isempty(match)
                                fragments(j) = all_deviants(j); % adds confirmed aberrant to list
                                aberrants = insertShape(double(aberrants),'circle',[match(1,1),match(1,2),3],...
                                    'color', 'green','LineWidth', 2); % confirms it on picture
                                
                                if(size(match,1) > 1) % if user clicked same chromosome multiple times, turns those green
                                    for k = 2:(size(match,1))
                                        aberrants = insertShape(double(aberrants),'circle',[match(k,1),match(k,2),3],...
                                            'color', 'green','LineWidth', 2); % confirms it on picture
                                    end
                                end
                                
                                % Only keeps track of user input that doesn't find a match
                                toRemove = find(ismember(user_input, match,'rows'));
                                if(~isempty(toRemove)) % in case user clicks on same chromosome twice
                                    user_input(toRemove, :) = [];
                                end
                            end
                        end
                        fragments(fragments == 0) = []; % remove all 0's
                        hold on, imshow(aberrants), hold off % show the user's clicks
                        
                        % invalid clicks are not allowed, restart if there are any
                        if(~isempty(user_input))
                            uiwait(msgbox('Not all user input points were found.','Warning','modal'));
                            close(gcf)
                        else
                            done = true;
                        end
                    end
                    
                    hold on, imshow(aberrants), hold off
                    pause(1) % allow user to see
                    close(gcf)
                    valid = true;
                    
                    % keep a list of fragments chunks, and redraw
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    found_fragments = unique([found_fragments, fragments]);
                    redraw = [redraw,true_aberrants(i)];
                case 'c' % chunk of a known broken fragment
                    chromosomes.PixelIdxList{true_aberrants(i)} = [];
                    chromosomes.NumObjects = chromosomes.NumObjects - 1;
                    valid = true;
                otherwise
                    uiwait(msgbox('Invalid Input.','modal'));
            end
        end
    end
    close(figure(1))
end
% clean up list
chromosomes.PixelIdxList = chromosomes.PixelIdxList(~cellfun('isempty',chromosomes.PixelIdxList))';

%% Spline Measuring - Redraw
% creates storage structure for data
image_data = struct('Chromosome_Length', {zeros(chromosomes.NumObjects,1)},...
    'Foci_Distances', {[]},'Chromosomes',{[]},'Male', isMale,...
    'Original', img, 'Composite', overlay, 'Dijkstra', overlay, 'Labeled', overlay);

for i = 1:length(redraw)
    done = false;
    while(~done) % redraw until satisfactory
        old_chromosome = logical(a); % creates blank logical matrix
        if(sum(found_fragments == redraw(i))) % if the redraw is a fragment, show all fragments
            for j = 1:length(found_fragments)
                old_chromosome(redFound.PixelIdxList{found_fragments(j)}) = 1;
            end
        else
            old_chromosome(redFound.PixelIdxList{redraw(i)}) = 1; % assigns all the listed pixels to 1
        end
        
        old_chromosome = imcomplement(old_chromosome);
        figure, imshow(old_chromosome)
        title('Please click along the length of one (1) of the chromosomes contained in the aberrant. e.g. Aberrants with 3 chromosomes within will be redrawn 3 times.')
        
        figh = gcf;
        if(i == 1), pos = get(figh,'position'); end
        set(figh,'position',resize_window(pos,6));
        
        % splining
        [y,x] = ginput;
        coord = [x,y];
        distance = 0;
        new_chromosome = logical(a); % creates blank logical matrix
        for j = 2:size(coord,1)
            % accumulate distances
            distance = distance + pdist([coord(j-1,:);coord(j,:)], 'euclidean');
            
            % create splines (1 per iteration)
            yLength = linspace(coord(j-1,2),coord(j,2),1000)';
            xLength = linspace(coord(j-1,1),coord(j,1),1000)';
            coords = unique(uint32(ceil([xLength,yLength])),'rows','stable');
            linear_indices = sub2ind(size(a),coords(:,1),coords(:,2));
            new_chromosome(linear_indices) = 1; % create the skinny line
        end
        
        % stores distance
        image_data.Chromosome_Length(length(chromosomes.PixelIdxList) + 1) = distance;
        
        % prepare to dilate
        [~, threshold] = edge(new_chromosome, 'sobel');
        new_chromosome = edge(new_chromosome,'sobel', threshold * fudgeFactor);
        
        % dilate and fill the drawn line
        seBegin = strel('line', 4, 100);
        seEnd = strel('line', 4, 0);
        lines = imdilate(new_chromosome, [seBegin seEnd]);
        new_chromosome = imfill(lines, 'holes');
        
        % smoothen
        newSmooth = imgaussfilt(double(new_chromosome), 1); % makes the image smooth so that there aren't cuts in the chromosome
        new_chromosome = imbinarize(newSmooth, graythresh(newSmooth)); % binarizes again
        
        % show
        show_new = logical(a);
        show_new(new_chromosome) = 1;
        
        imshowpair(old_chromosome,imcomplement(show_new))
        figh = gcf;
        set(figh,'position',resize_window(pos,6));
        
        choice = questdlg('Is this redrawn chromosome okay?','Redraw Confirmation','Yes','Redraw','Yes');
        if(strcmp(choice,'Yes'))
            done = true;
        end
    end
    
    % store
    chromosomes.PixelIdxList{length(chromosomes.PixelIdxList) + 1} = find(new_chromosome ~= 0);
    close(gcf)% store the new chromosome
end
close(gcf)

%% Measurements?
% pre-allocate storage for each cell's outline to be used later
outlines = cell(1, chromosomes.NumObjects);
% remembers the chromosomes we already have measures for
redrawn = chromosomes.NumObjects - (length(redraw)-1):chromosomes.NumObjects;
for i = 1:chromosomes.NumObjects % isolates each chromosome found
    
    blank = logical(a); % creates blank logical matrix
    blank(chromosomes.PixelIdxList{i}) = 1; % assigns all the listed pixels to 1
    
    % outline
    [~, threshold] = edge(blank, 'sobel');
    perim = edge(blank,'sobel', threshold * fudgeFactor);
    outlines(i) = {double(perim)}; % save outline for later use
    
    if(~sum(i == redrawn)) % if it's not an aberrant, store the length
        perimSum = sum(sum(perim))/2;
        image_data.Chromosome_Length(i) = perimSum - 5; % may need tweaking
    end
    % subtracting five is chosen arbitrarily as a correction for the slight
    % curve towards the tip of a chromosome
end
if(exist('XY','var')) % saves XY outline if it needs to
    blank = logical(a); % creates blank logical matrix
    blank(XY.PixelIdxList) = 1;
    % outline
    [~, threshold] = edge(blank, 'sobel');
    perim = edge(blank,'sobel', threshold * fudgeFactor);
    outlines(chromosomes.NumObjects + 1) = {double(perim)}; % save outline for later use
    chromosomes.NumObjects = chromosomes.NumObjects + 1;
    chromosomes.PixelIdxList{chromosomes.NumObjects} = XY.PixelIdxList;
end
%% Dijkstra's Approach to automated measuring of inter-foci distance
% create a map of the image to pass to Dijkstra's method
overall = overlay; % copies the type from a pre-existing variable
overall(:) = 0; % clears it out
no_centromere = zeros(chromosomes.NumObjects,1);
chromosome_centers = cell(chromosomes.NumObjects,1);
for i = 1:chromosomes.NumObjects
    body = logical(a); % blank map of 0's
    body(chromosomes.PixelIdxList{i}) = 1; % body of chromosome
    map = 2*double(body); % begins map creation
    skeleton = double(bwmorph(map,'skel',Inf)); % saves skeletal chromosome
    
    map(outlines{i} ~= 0) = 3; % makes the perimeter higher cost
    map(skeleton == 1) = 1; % makes the skeletal path the lowest cost
    
    % create body of chromosome for Dijkstra image
    blank = logical(a);
    blank(chromosomes.PixelIdxList{i}) = 1;
    
    % calculate center so we can later label it
    tmp = struct2cell(regionprops(blank, 'centroid'));
    chromosome_centers{i} = cellfun(@(x) uint32(x),tmp,'UniformOutput',false);
    
    % Dijkstra iamge
    overall = overall + cat(3, blank, logical(skeleton), logical(outlines{i}));
    %     imshowpair(img, overall, 'montage'); hold on
    
    % store Dijkstra image
    image_data.Dijkstra = logical(overall);
    
    % find centroid of centromere
    centromere = regionprops(body & bwB, 'centroid');
    centromere = struct2cell(centromere);
    % won't look for centromeres in an XY chromosome
    
    if(~(exist('XY','var') && i == chromosomes.NumObjects))
        % if an error occurs here, the centromere was missed when redrawing
        try
            centromere = uint32(centromere{1}); % unpack the structure
        catch
            warning('A non-XY chromosome (#%i) with foci present could not be measured due to lack of centromere.',i)
            no_centromere(i) = i;
            continue % :(
        end
    end
    
    % find centroids of foci
    foci = regionprops(body & bwG, 'centroid');
    if(isempty(foci)), continue, end % skip this chromosome if there aren't any foci
    foci = int32(cat(1, foci.Centroid));
    
    
    if(exist('XY','var') && i == chromosomes.NumObjects)
        [distance, ~] = dijkstra_image(map, XY.Centromeres(1,:), XY.Centromeres(2,:));
        XY.Length = distance;
        image_data.Chromosome_Length(i) = distance;
        image_data.Foci_Distances{i} = -1;
        hold off
        continue
    end
    [distance, ~] = dijkstra_image(map, centromere, foci);
    image_data.Foci_Distances{i} = sort(distance);
    hold off
end

image_data.Chromosomes = chromosomes.PixelIdxList;
if(size(image_data.Foci_Distances, 1) == 1)
    image_data.Foci_Distances = image_data.Foci_Distances'; % makes data organization consistent
end
close all

% display for user
newFound = struct();
newFound.('Connectivity') = 8;
newFound.('ImageSize') = size(a);
newFound.('NumObjects') = length(image_data.Chromosome_Length);
newFound.('PixelIdxList') = image_data.Chromosomes';

subplot(2,2,1), imshow(img), title('Original Image')
subplot(2,2,2), imshow(double(overlay)), title('Composite Image')

% label the chromosomes with their number
color_map = label2rgb(labelmatrix(newFound),'hsv','w', 'shuffle');
for i = 1:length(chromosome_centers)
    center = chromosome_centers{i};
    color_map = insertText(color_map,center{1}, i);
end
image_data.Labeled = color_map; % store
subplot(2,2,3), imshow(color_map), title('Labeled Color Map')
subplot(2,2,4), imshow(overall), title('Final Dijkstra Image')

figh = gcf;
pos = get(figh,'position');
set(figh,'position',resize_window(pos,8));

%% Extract Data -> [Chromosome Length, Foci Distances (if any)]
new_size = max(cellfun('length',image_data.Foci_Distances));
extracted_data = zeros(length(image_data.Foci_Distances),new_size);
for i = 1:length(image_data.Foci_Distances)
    new_row = zeros(1,new_size);
    for j = 1:length(cell2mat(image_data.Foci_Distances(i)))
        values = cell2mat(image_data.Foci_Distances(i));
        new_row(j) = values(j);
    end
    extracted_data(i,:) = new_row;
end
sort(extracted_data, 2); % arrange foci (column wise) by distance

% create final table
extracted_data = horzcat(image_data.Chromosome_Length,extracted_data,-no_centromere);
names = horzcat('Length', strcat(repmat({'Foci_'},size(extracted_data,2)-2,1),...
    cellfun(@(x) num2str(x),num2cell(1:size(extracted_data,2)-2)'))','Skipped');
extracted_data = array2table(extracted_data, 'VariableNames', names);
% display(extracted_data)
end