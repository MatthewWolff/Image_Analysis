%% Image Analysis Batch Runner
location = '~/Desktop/College/Research/PayseurLab/';
files = dir(strcat(location,'/*.tif')); % finds all matching files
while(isempty(files)) % user input when necessary
    location = inputdlg('Please supply the filepath to the image folder');
    location = location{1};
    files = dir(strcat(location,'/*.tif')); % finds all matching files
end
file_names = {files.name}';
raw_names = cellfun(@(x) x(1:(length(x)-4)),file_names,'UniformOutput', false); % 4 == .tif length
file_location = {location};
file_paths = strcat(file_location,file_names);

% easily create table and then convert to data structure
extracted_data = zeros(length(file_paths),1);
image_data = zeros(length(file_paths),1);
batch_data = table(file_names,extracted_data,image_data,'RowNames',raw_names);
batch_data = table2struct(batch_data);

for i = 1:length(file_paths) % this should go so that when something goes wrong, it continues
    
    try
        [extracted,image] = image_analysis(file_paths{i});
        batch_data(i).extracted_data = extracted;
        batch_data(i).image_data = image;
    catch ERROR
        error_report = getReport(ERROR);
        continue
    end
end
