function MSK2NPYsingle(NPYfolder,freeformmask)

filename_MSK = freeformmask.filename;
name = freeformmask.name;
%% Set intial parameter
main_trans = 0; % main feature transmittance 
main_phas = 0; % main feature phase
main_group = 0; % main feature group(0)
main_poly_info = [main_trans main_phas main_group];
% pixel_size_mask = 32; % nano meter
% pixel_size_mask = 40; % nano meter
pixel_size_mask = 5; % nano meter
Trainsize = 896; % Training image size
%xTrainImages = cell(1,length(filename_MSK));

% xTrainImages = cell(1,2000);
% XTrain = ones(Trainsize,Trainsize,1,2000);


%% Read mask
% target_filename_MSK = fullfile(path,filename_MSK{1,i});
[~,sim_region,back_ground_info,msk_poly_target,msk_poly_info,n_poly_target] = Parse_MSKfile(filename_MSK); % parse MSK information
% [Main_Poly,Main_Poly_Info,~] = Main_Polygon(msk_poly_target,msk_poly_info,n_poly_target,sim_region);
[sraf_poly_target] = extract_sraf(msk_poly_target,msk_poly_info,n_poly_target,sim_region);
n_poly_target = length(sraf_poly_target);
%% Pixel based mask
[target_pixelbased,~,~,~,~,~,~] = Find_Mask_Group_Center_Pixel(back_ground_info,sim_region,pixel_size_mask,sraf_poly_target,n_poly_target,main_poly_info);% polygon to pixel-based
% target_pixelbased = target_pixelbased(1:128,1:128);
% imagesc(abs(target_pixelbased));
XTrain = ones(Trainsize,Trainsize);

% Padding
if mod(size(target_pixelbased,1),2) == 0
    top_pad = (Trainsize - size(target_pixelbased,1))/2;
    bottom_pad = (Trainsize - size(target_pixelbased,1))/2;
else
    top_pad = floor((Trainsize - size(target_pixelbased,1))/2);
    bottom_pad = floor((Trainsize - size(target_pixelbased,1))/2)+1;
end

if mod(size(target_pixelbased,2),2) == 0
    left_pad = (Trainsize - size(target_pixelbased,2))/2;
    right_pad = (Trainsize - size(target_pixelbased,2))/2;
else
    left_pad = floor((Trainsize - size(target_pixelbased,2))/2);
    right_pad = ceil((Trainsize - size(target_pixelbased,2))/2);
end
XTrain(:,:) = back_ground_info(1,1) * XTrain(:,:);
XTrain(top_pad+1:(Trainsize-bottom_pad),left_pad+1:(Trainsize-right_pad)) ...
= abs(target_pixelbased(:,:));

imagedata(:,:) = imbinarize(XTrain(:,:));
imagedata = double(imagedata);
% name = regexp(filename_MSK{1,1},'_','split');
% name = regexp(filename_MSK,'.M','split');
save_name = [name '.npy'];
writeNPY(imagedata, fullfile(NPYfolder,save_name));
end