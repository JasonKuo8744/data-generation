%% load all data (source and target file) (constant parameter)
close all;
clear ;
clc;
%%
folderworkspace = pwd;

% Choose preparation folder
PathName = uigetdir(folderworkspace,'Specify Prepare Folder');
if PathName == 0
    return;
end

% Choose result folder 
folderdestination = uigetdir(folderworkspace,'Specify Results Folder');
if folderdestination == 0
  return;
end

% Choose NPY folder 
SRAFpathforNPY = uigetdir('.npy','Specify SRAF NPY Folder');
if SRAFpathforNPY == 0
    return;
end

% Choose NPY folder 
EPEpathforNPY = uigetdir('.npy','Specify EPE NPY Folder');
if EPEpathforNPY == 0
    return;
end

% Choose NPY folder 
PrintOutPath = uigetdir('.npy','Specify Print Out Folder');
if PrintOutPath == 0
    return;
end

listing = dir(PathName);
name = {listing(1:end).name};

index = strfind(name,'.pl2');
for i = 1:length(index)
    flag = isempty(index{1,i});
    if flag == 0
        number = i;
    end
end
filename_pl2{1,1} = name{1,number};


index1 = strfind(name,'.mat');
for i = 1:length(index1)
    flag = isempty(index1{1,i});
    if flag == 0
        number = i;
    end
end
filename_mat{1,1} = name{1,number};

index2 = strfind(name,{'.xlsx'});
for i = 1:length(index2)
    flag = isempty(index2{1,i});
    if flag == 0
        number = i;
    end
end
index2 = name{1,number};
index2 = strsplit(index2,'$');
filename_xlsx{1,1} = index2{1,2};

index3 = strfind(name,'.SRC');
for i = 1:length(index3)
    flag = isempty(index3{1,i});
    if flag == 0
        number = i;
    end
end
filename_SRC{1,1} = name{1,number};

index4 = strfind(name,'.MSK');
for i = 1:length(index4)
    flag = isempty(index4{1,i});
    if flag == 0
        number = i;
    end
end
filename_MSK{1,1} = name{1,number};

% filename_pl2 = name(contains(name,'.pl2')==1);
% filename_mat = name(contains(name,'.mat')==1);
% filename_xlsx = name(contains(name,{'.xlsx','.xls'})==1);
% filename_SRC = name(contains(name,'.SRC')==1);
% filename_MSK = name(contains(name,'.MSK')==1);

if length(filename_mat)~=1
    filename_mat = filename_mat(contains(filename_mat,'cs','IgnoreCase',true));
end
% filename_xlsx = filename_xlsx(~contains(filename_xlsx,'~$'));

if length(filename_pl2) == 1 && length(filename_mat) == 1 && length(filename_xlsx) == 1 && length(filename_SRC) == 1 && length(filename_MSK) == 1
    
    filename_clip = fullfile(PathName,filename_pl2);
    filename_CS = fullfile(PathName,filename_mat);
    filename_setting = fullfile(PathName,filename_xlsx);
    source_filename_SRC = fullfile(PathName,filename_SRC);
    target_filename_MSK = fullfile(PathName,filename_MSK);
    
else
    return
end

clipsData = extractclipdata(filename_clip,filename_CS);
if isempty(clipsData)
  return;
end

settingData = extractsettingdata(filename_setting{1,1},clipsData.n_clip);
lengthdata.n_EPE_fitness = 9;
lengthdata.n_variationcondition = 1;
settingData.ParamProcessCondition.lengthdata = lengthdata;
% cssamplelength = 40;
% cssamplesiz = 1;
% clipsData = addCSline(clipsData,cssamplelength,cssamplesiz);

freeformdata.name = filename_SRC{1,1}; % SRC nameen
freeformdata.filename = source_filename_SRC{1,1};
%% in extractsourcetargetdata change the transmittance and phase of main feature and background  
t0 = clock;
lithoData = extractsourcetargetdata(source_filename_SRC{1,1},target_filename_MSK{1,1},settingData.ParamSimulationCondition);
TargetData = lithoData.TargetData;
settingData.ParamProcessCondition.flag_pOPC = TargetData.flag_pOPC;
settingData.ParamProcessCondition.OPCData = lithoData.OPCData;

% figure;
% pupilW = lithoData.IlluminateData.pupilW;
% surf(pupilW);

%% clipdata established
% tic %timer1

cssamplelength = 40;
segmentsize = 100;
lineEndLength = 50;
jog = 0.5;
damping = 0.7;
CSdata = add_csdata(TargetData,segmentsize,cssamplelength,lineEndLength,jog);
clipsData.CSdata = CSdata;
clipsData.filename_clip = filename_clip;
[~, filenamesave_clip, ~] = fileparts(filename_clip{1});
filenamesave_clip = {filenamesave_clip};
clipsData.filenamesave_clip = filenamesave_clip;
clipsData.main_poly = TargetData.main_poly;
clipsData.n_main = length(TargetData.main_poly);
clipsData.main_center = TargetData.main_center;
clipsData.n_clip = length(filename_clip);

%% kmeans block
% Point_Segment = 10; % Segment size
% Point_Interval = 20; % Interval size
% [Polygon_Boundary_Point,~] = Find_Polygon_Point(length(TargetData.main_poly),TargetData.main_poly,TargetData.sim_region,Point_Segment,Point_Interval); % find the point on the polygon edge
% Polygon_Boundary_Point.bound = Polygon_Boundary_Point.bound';
% Polygon_Boundary_Point.poly = Polygon_Boundary_Point.poly';

%% calculate interference map
source = lithoData.SourceData.source;
NA     = lithoData.IlluminateData.NA;
lambda = lithoData.IlluminateData.lambda;
pupil  = lithoData.IlluminateData.pupil;
pupil_de = lithoData.IlluminateData.pupil_de;
pupil_mde = lithoData.IlluminateData.pupil_mde;
target_pixelbased = lithoData.TargetData.target_pixelbased;
pixel_size_mask = TargetData.pixel_size_mask;

[IMap.original] = litho2DTCC(source,pupil,NA,lambda,target_pixelbased,pixel_size_mask);
% [IMap1.original] = litho2DTCC(source,pupil_de,NA,lambda,target_pixelbased,pixel_size_mask);
% [IMap2.original] = litho2DTCC(source,pupil_mde,NA,lambda,target_pixelbased,pixel_size_mask);
[IMap.normalize,IMap.abs,IMap.phase,TargetData.flag_vertical] = NormalizeIM(IMap.original,TargetData,IMap.original); % let interference map value (0~1)
% [IMap1.normalize,IMap1.abs,IMap1.phase,TargetData.flag_vertical] = NormalizeIM(IMap1.original,TargetData,IMap1.original);
% [IMap2.normalize,IMap2.abs,IMap2.phase,TargetData.flag_vertical] = NormalizeIM(IMap2.original,TargetData,IMap2.original);

%% Draw IM map
% figure;
% xmax = TargetData.mask_dimension(2);
% xmin = TargetData.mask_dimension(4);
% ymax = TargetData.mask_dimension(1);
% ymin = TargetData.mask_dimension(3);
% x = xmin:TargetData.pixel_size_mask:xmax;
% y = ymin:TargetData.pixel_size_mask:ymax;
% imagesc(x,y,abs(IMap.normalize));
% colormap('jet');
% colorbar;
% axis image;
% axis xy;
% set(gca,'FontName', 'Times New Roman')
% set(gca,'FontSize', 18)
% set(gca,'FontWeight','Bold')
% xlabel('X(nm)');
% ylabel('Y(nm)');

%% collect SRAF rewrite the mask in prolith and initial genetic algorithm parameter (user decide)
flag_main_positive = TargetData.flag_main_positive;
if flag_main_positive
    Precheck_IMmap.level_min = 0; % ratio to get SRAF (x_level = y_level) #
    Precheck_IMmap.level_max = 0.5; % #
    level_min_map = +(IMap.normalize > Precheck_IMmap.level_min);
    level_max_map = +(IMap.normalize > Precheck_IMmap.level_max);
else
    Precheck_IMmap.level_min = 0.5; % ratio to get SRAF (x_level = y_level) #
    Precheck_IMmap.level_max = 1; % #
    level_min_map = +(IMap.normalize < Precheck_IMmap.level_min) | +(abs(IMap.normalize-Precheck_IMmap.level_min) < 0.00001);
    level_max_map = +(IMap.normalize < Precheck_IMmap.level_max) | +(abs(IMap.normalize-Precheck_IMmap.level_min) < 0.00001);
end
  
% figure(1);
% imagesc(level_min_map);
% axis image;
% axis xy;
% figure(2);
% imagesc(level_max_map);
% axis image;
% axis xy;

% max_matching_dis = 3*TargetData.pixel_size_mask; % identify the characters
% max_move_dis = 7*TargetData.pixel_size_mask; % maximum distance of moving polygon
%%
max_section_box = 5*TargetData.pixel_size_mask; 
[section_map] = dividemaptosection(TargetData,max_section_box);
clipsData.main_poly = section_map.main_poly;
clipsData.n_main = section_map.n_main;
clipsData.main_center = section_map.main_center;
clipsData.TargetData = TargetData;

%% Draw section map
% figure;
% drawPolygon(section_map.poly_me)
% % drawPolygon(section_map.section_sraf_poly,'g')
% drawPolygon(section_map.main_poly,'b')
% drawPolygon(section_map.poly_va_ver,'g','lineWidth', 2)
% axis image
% set(gca,'Fontsize',20)
% xlabel('nm','Fontsize',20)
% ylabel('nm','Fontsize',20)
%%
% GAP.nVar = TargetData.target_pixelnumber(1,1)*TargetData.target_pixelnumber(1,2); % number of decision variables
% GAP.VarSize = [1 GAP.nVar]; % decision variables matrix size
GAP.nPop = 20; % initial population size #
k_number = 0; % initial IM SRAF
n_random = 20; % initial random SRAF
GAP.first_ten_level = Precheck_IMmap.level_min:((Precheck_IMmap.level_max-Precheck_IMmap.level_min)/k_number):Precheck_IMmap.level_max;

%% prepare 
Mask_PixelMap = cell(k_number,1);

empty_individual.Polygon_Vertex = [];
empty_individual.Polygon_Information = []; % [xmin xmax ymin ymax sraf/main trans(background:0/feature:1)]
empty_individual.Feature_Group = [];
empty_individual.GA_Gene = []; % [sraf/main section_number gene]

empty_individual.Cost = [];
empty_individual.EPE = [];

empty_individual.Main_Number_Right = [];

First_Pop = repmat(empty_individual,k_number,1);
Random_Pop = repmat(empty_individual,n_random,1);
Random_F_Pop = repmat(empty_individual,n_random,1);
pop = repmat(empty_individual,k_number,1);

unchange_poly = section_map.unchange_poly;
unchange_poly_info = section_map.unchange_poly_info;
main_poly = section_map.main_poly;
main_poly_info = section_map.main_poly_info;

%% connect PROLITH
[PROLITHdata] = initializeprolithPPI(clipsData,settingData.ParamSimulationCondition);
if isempty(PROLITHdata)
  return;
end
PROLITHdata.freeformdata = freeformdata;
%% i.c. random sraf + target
flag_main_positive = TargetData.flag_main_positive;
random_Sraf = section_map.gene_sraf_group;
[~,sorts] = sort([random_Sraf.n_gene],'descend');
random_Sraf = random_Sraf(sorts);
n_sraf_all = sum([random_Sraf.n_gene]);
n_sraf_random = randi(n_sraf_all,n_random,1);
for count = 1:n_random
%  for count = 1:n_random
     ind_sraf_random = randsample(n_sraf_all,n_sraf_random(count));
     poly_rinfo = zeros(n_sraf_random(count),6);
     GA_Gene = cell(n_sraf_random(count),1);
     
     for kcount = 1:n_sraf_random(count)
         ind = ind_sraf_random(kcount);
         GA_Gene{kcount} = [1 sorts(ind) 1];
         info = random_Sraf(ind).gene{1};
         poly_rinfo(kcount,:) = info;
     end
     
     [index1,index2] = meshgrid(1:n_sraf_random(count)-1,1:n_sraf_random(count));
     index = [reshape(index1,[],1) reshape(index2,[],1)];
     index(index(:,1)>=index(:,2),:)=[];
     check_poly = [];
     
     for number = 1:length(index(:,1))
         
         index1 = index(number,1);
         x1min = poly_rinfo(index1,1) - 16.25/2 ;
         x1max = poly_rinfo(index1,2) + 16.25/2 ;
         y1min = poly_rinfo(index1,3) - 16.25/2 ;
         y1max = poly_rinfo(index1,4) + 16.25/2 ;
         
         index2 = index(number,2);
         x2min = poly_rinfo(index2,1) - 16.25/2 ;
         x2max = poly_rinfo(index2,2) + 16.25/2 ;
         y2min = poly_rinfo(index2,3) - 16.25/2 ;
         y2max = poly_rinfo(index2,4) + 16.25/2 ;
         
         if ~(x1max < x2min || x1min > x2max || y1max < y2min || y1min > y2max)
             if max(abs(x1max - x1min),abs(y1max - y1min)) > max(abs(x2max - x2min),abs(y2max - y2min))
                 check_poly = [check_poly ; index1 index2];
             else
                 check_poly = [check_poly ; index2 index1];
             end
         end
     end
     
     index1 = [];
     index2 = [];
     indexb = [];
     n=0;
     check_poly1 = check_poly;
     while ~isempty(check_poly1)
         
         n = n + 1;
         check_siz = length(check_poly1);
         kcount = find(unique(check_poly1));
         if sum(histcounts(check_poly1,1:n_sraf_random(count)+1)) ~= kcount(end)
             [~,index1(n,1)] = max(histcounts(check_poly1,1:n_sraf_random(count)+1));
             a = find(sum(check_poly1==index1(n,1),2)==1);
             indexb = [indexb ; check_poly1(a,2)];
             index2 = mod(find(check_poly1 == index1(n,1)),check_siz);
             index2(index2==0) = check_siz;
             check_poly1(index2,:) = [];
         else
             if flag_main_positive == 2
                 indexb = [indexb;check_poly1(:,2)];
                 check_poly1 = [];
             else
                 index1 = [index1;check_poly1(:,2)];
                 check_poly1 = [];
             end
             break
         end
     end
     
     if flag_main_positive == 2
         if ~isempty(indexb)
             GA_Gene(indexb) = [];
             poly_rinfo(indexb,:) = [];
         end
     else
         if ~isempty(index1)
             GA_Gene(index1) = [];
             poly_rinfo(index1,:) = [];
         end
     end
     
     Random_Pop(count).GA_Gene = GA_Gene;
     Random_Pop(count).Polygon_Information = poly_rinfo;
     name_first = ['Random_Pop' num2str(count)];
     
     [Random_poly,Random_info] = polyinfo2vertex(Random_Pop(count),TargetData);
     poly = [main_poly;
         unchange_poly;
         Random_poly];
     info = [main_poly_info;
         unchange_poly_info;
         Random_info];
     Random_Pop(count).Name = name_first;
     % >> Do own OPC
     if ~settingData.ParamProcessCondition.flag_pOPC
         if TargetData.flag_1D == 0
            [OPC_mask,OPC_pop] = model_based_OPC(poly,info,clipsData,PROLITHdata,settingData,folderdestination,'iteration',20,'damping factor',damping);
         else
             [OPC_mask,OPC_pop] = model_based_1DOPC(poly,info,clipsData,PROLITHdata,settingData,folderdestination,'iteration',20,'damping factor',damping);
         end
         Random_Pop(count).EPE = OPC_mask.EPE;
         Random_Pop(count).Cost = OPC_mask.Cost;
     else
         freeformmask = writemskfile(Random_Pop(count),TargetData,section_map,name_first,folderdestination);
         PROLITHdata.freeformmask = freeformmask;
         clipsData.filenamesave_clip{1} = freeformmask.name;
         settingData.ParamProcessCondition.flag_save_opc = false;
         settingData.ParamProcessCondition.iteration = 0;
         [Random_Pop(count).EPE] = evaluateclips(clipsData,PROLITHdata,settingData.ParamProcessCondition,folderdestination);
         [Random_Pop(count).Cost] = calculatefitnessvalue(Random_Pop(count).EPE,true,settingData.ParamWeightCondition.weighting_EPE);
         EPE_name_first = ['EPE_Random_Pop' num2str(count) '.npy'];
         Check_name_first = ['Check_Random_Pop' num2str(count) '.npy'];
%          delete([fullfile(folderdestination,name_first) '.pl2']);
         MSK2NPYsingle(SRAFpathforNPY,freeformmask);
%          writeNPY(Random_Pop(count).Cost, fullfile(EPEpathforNPY,EPE_name_first));
         
%          if sum(Random_Pop(count).EPE.sraf_record ~= 0) == 0
%              print_out = 1; % SRAF print out
%              writeNPY(print_out, fullfile(PrintOutPath,Check_name_first));
%          else
%              print_out = 0; % SRAF didn't print out
%              writeNPY(print_out, fullfile(PrintOutPath,Check_name_first));        
%          end
% 
     end
%      disp(['Random Population ' num2str(count) ' Cost ' num2str(Random_Pop(count).Cost)]);
%      Random_Pop = repmat(empty_individual,n_random,1);
 end
 
 
 
 
 %% i.c. generated by imap 
for number = 1:k_number
    if flag_main_positive
        Mask_PixelMap{number} = IMap.normalize > GAP.first_ten_level(number); % map is logical (maintran > backtran)
    else
        Mask_PixelMap{number} = IMap.normalize < GAP.first_ten_level(number); % map is logical (maintran < backtran)
    end
    % Calculate SRAF position, reshape into rectangular and assign gene
    [First_Pop(number),section_map] = createfirstpopulation(Mask_PixelMap{number},First_Pop(number),TargetData,section_map,16.25); 
    main_number_right = First_Pop(number).Main_Number_Right;

    if main_number_right
        
        [IMAP_poly,IMAP_info] = polyinfo2vertex(First_Pop(number),TargetData);
        poly = [main_poly;
            unchange_poly;
            IMAP_poly];
        info = [main_poly_info;
            unchange_poly_info;
            IMAP_info];
        
        name_first = ['First_Pop' num2str(number)];
        First_Pop(number).Name = name_first;
        if ~settingData.ParamProcessCondition.flag_pOPC
            if TargetData.flag_1D == 0
                [OPC_mask,OPC_pop] = model_based_OPC(poly,info,clipsData,PROLITHdata,settingData,folderdestination,'iteration',20,'damping factor',0.6);
            else
                [OPC_mask,OPC_pop] = model_based_1DOPC(poly,info,clipsData,PROLITHdata,settingData,folderdestination,'iteration',20,'damping factor',0.6);
            end
            First_Pop(number).EPE = OPC_mask.EPE;
            First_Pop(number).Cost = OPC_mask.Cost;
        else
            freeformmask = writemskfile(First_Pop(number),TargetData,section_map,name_first,folderdestination);
            PROLITHdata.freeformmask = freeformmask;
            clipsData.filenamesave_clip{1} = freeformmask.name;
            settingData.ParamProcessCondition.flag_save_opc = false;
            settingData.ParamProcessCondition.iteration = 0;
            [First_Pop(number).EPE] = evaluateclips(clipsData,PROLITHdata,settingData.ParamProcessCondition,folderdestination);
            [First_Pop(number).Cost] = calculatefitnessvalue(First_Pop(number).EPE,true,settingData.ParamWeightCondition.weighting_EPE);
            EPE_name_first = ['EPE_First_Pop' num2str(number) '.npy'];
            Check_name_first = ['Check_First_Pop' num2str(number) '.npy'];
%             delete([fullfile(folderdestination,name_first) '.pl2']);
            MSK2NPYsingle(SRAFpathforNPY,freeformmask);
            writeNPY(First_Pop(number).Cost, fullfile(EPEpathforNPY,EPE_name_first));
            
            if sum(First_Pop(number).EPE.sraf_record ~= 0) == 0
                print_out = 1; % SRAF didn't print out
                writeNPY(print_out, fullfile(PrintOutPath,Check_name_first));
            else
                print_out = 0; % SRAF print out
                writeNPY(print_out, fullfile(PrintOutPath,Check_name_first));        
            end
            
        end
    else
        First_Pop(number).Cost = Inf;
    end
    disp(['Imap Population ' num2str(number) ' Cost: ' num2str(First_Pop(number).Cost)]);
    First_Pop = repmat(empty_individual,k_number,1);
end

%%
% MSK2NPY(folderdestination,SRAFpathforNPY);

