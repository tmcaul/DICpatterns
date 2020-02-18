clear all

direc='/Users/tom/Documents/GitHub/tpm_DICpatterns/Patterns_Tianhong';
%dir='D:\Tom\Box Sync\Patterns\Patterns_120319'

mkdir(direc)
cd(direc)

addpath(genpath('/Users/tom/Documents/GitHub/tpm_DICpatterns'));

% Good idea to stitch together multiple images / fields of view.

% Resolution of DIC is limited by subset size. Need a no. of speckles per
% subset (ideally 3x3), and need speckles to be large enough (3x3) pixels in the images
% that are taken.

% Therefore, stitching images allows finer speckles, increased density, more/smaller subsets, and
% better resolution.

% Define a speckle size to recommend a HFW for stitch images. This is
% because size of speckles needs to be an independent variable for
% stability and contrast purposes.

% REMINDER:
% DENSITY of speckles determines subset size AND THEREFORE RESOLUTION
% SIZE of speckles determines stitch image size (and therefore how many needed to cover a reigon, as speckle should take up 3x3 pixels)

region_x=5; % dimensions of region to be patterned in um
region_y=10;

grids=0; %grids on or off
grid_spacing=10;%grid dimension in um
grid_width=0.1; %;%grid

speckles=1; %speckles on or off
speckle_width=0.147; %in um
subset_dimension=0.7; %square subset dimension in um: There should be 3x3 speckles per subset. ie 9 per this no squared.
req_speckles=5; %no of required speckles per subset: 9? 5?

% speckle pattern generation parameters:
no=3; %no of speckles that affect shifts
%too high and av_distance gets too big; too low and not enough shift.
thresh_high=0.4; %of +/- of equilibrium separation that will shift
thresh_low=0.4;
mag=0.2; %of shift of micro-adjustments
its=8000; %no of iterations for micro-adjustments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%speckle shape for output image
shape=...
[0,0,1,0,0;
0,1,1,1,0;
1,1,1,1,1;
0,1,1,1,0;
0,0,1,0,0;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stitch image parameters:
img_size=[1768, 2048]; %rows and columns so x, y reversed [884, 1024]
d_um2=req_speckles./(subset_dimension^2); %density of speckles: no. per um^2
pixel_width=speckle_width./3; %3 pixels per speckle determines recommended image size (in um per pixel)
hfw=img_size(2).*pixel_width; % j dimension is x - determines recommended hfw of stitch images.
%
% Initialise things in terms of pixels rather than space
region_i=round(region_y./pixel_width);
region_j=round(region_x./pixel_width);
d_pixels2=d_um2.*(pixel_width.^2); %density of speckles in no. per sq. pixel.= no. per um2 * pixelsperum^2
speckle_separation=round(sqrt(1/d_pixels2)); %speckle separation in pixels

%%%%%%%%%%%%%%%%%%%
% Set up region
%%% set up lattice of speckle centres with density set by d (and therefore speckle separation)
%%% each speckle then can move by a random number of pixels in x and y
%%% amount by which each pixel can move is 1/3 of ordered distance

% region is pattern without shape convolution
region=zeros(region_i,region_j);

if speckles==1;
for i=1:region_i;
    for j=1:region_j;
        
        % random offset between +/- 1/3 of speckle separation
        if mod(i,speckle_separation)==0 & mod(j,speckle_separation)==0;
            
           offset1=randi([-round(speckle_separation/3), round(speckle_separation/3)]);
           offset2=randi([-round(speckle_separation/3), round(speckle_separation/3)]);
           region(i+offset1-round(speckle_separation/3),j+offset2-round(speckle_separation/3))=1; 
        end
        j=j+1;
    end
    i=i+1;
end

clear offset1 offset2 i j a speckles_init 
original_region=region;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimise speckle pattern
%%%% do micro-adjustments on random speckles
f=waitbar(1./its, 'optimising...');

for index=1:its;
    waitbar(index./its);

    %[D,idx]=bwdist(image);
    [speckle_i,speckle_j]=find(region); %find non-zero elements and their indices
    speckle_locs=[speckle_i, speckle_j];
    clear speckle_i speckle_j
    clear distances distances_sorted distances_sorted_cropped speckles_closest sq_distances vectors net_vector loc speckle_pos_new;
    
% trial a random pixel and get closest speckles
    trial_i=randi([30,region_i-30]);
    trial_j=randi([30,region_j-30]);
    loc=[trial_i,trial_j];
    clear trial_i trial_j
    
    distances=speckle_locs-loc; %distances between trial location and closest speckles
    sq_distances=sqrt((distances(:,1)).^2+(distances(:,2)).^2); %square distances
    [distances_sorted(:,1),distances_sorted(:,2)]=sort(sq_distances); %outputs min. square distances and the index in the array
    distances_sorted_cropped=distances_sorted(1:no,:); %crops array to the speckle itself and no. closest
    
    for b=1:no; %generate closest speckles i,j indices
        speckles_closest(b,:)=speckle_locs(distances_sorted(b,2),:);
    end
    
    av_distance=mean(distances_sorted_cropped(:,1));

    %%%%% make adjustments to speckles
    if av_distance>(1+thresh_high)*speckle_separation;
        %move closest speckles to be slightly closer to trial location
        for b=1:no; %consider the requisite number of speckles specified at top
            vector=loc-speckles_closest(b,:); %vector from closest speckle to trial location
            region(speckles_closest(b,1),speckles_closest(b,2))=0; %remove considered speckle
            region(speckles_closest(b,1)+round(mag.*vector(1)),speckles_closest(b,2)+round(mag.*vector(2)))=1; %move considered speckle a bit closer
        end
    end
    
    if av_distance<(1-thresh_low)*speckle_separation;
        %move closest speckles to be slightly further from trial location
        for b=1:no; %consider the requisite number of speckles specified at top
            vector=loc-speckles_closest(b,:); %vector from closest speckle to trial location
            if speckles_closest(b,1)-round(mag.*vector(1))>0&speckles_closest(b,2)-round(mag.*vector(2))>0;
                region(speckles_closest(b,1),speckles_closest(b,2))=0; %remove considered speckle
                region(speckles_closest(b,1)-round(mag.*vector(1)),speckles_closest(b,2)-round(mag.*vector(2)))=1; %move considered speckle a bit further
            end
        end
    end
    
end
close(f)

clear a b speckles_closest index speckle_locs speckle_separation sq_distances loc speckle_pos_new c distances distances_sorted distances_sorted_cropped net_vector vectors f vector av_distance

%%%%%%%%%%%%%%%%
% convolution of shape speckle at each centre (3x3 pixels optimal)
pattern=conv2(shape,region);
pattern_1=conv2(shape,original_region); % in case want to compare optimised / unoptimised

% Final speckle locations:
% (Remember i and j are inverted x and y...)
[final_speckle_i,final_speckle_j]=find(region);
speckle_x=final_speckle_j.*pixel_width;
speckle_y=final_speckle_i.*pixel_width;

else
    speckle_x=0;
    speckle_y=0;
end

% add grids
if grids==1;
    grid_separation_pixels=round(grid_spacing./pixel_width);
    grid_width_pixels=round(grid_width./pixel_width);

    for i=1:region_i;
        for j=1:region_j;
            if mod(i,grid_separation_pixels)==0||mod(j,grid_separation_pixels)==0;   
                pattern(i+grid_width_pixels,j+grid_width_pixels)=1;  
            end
            j=j+1;
        end
        i=i+1;
    end
end

clear i j grid_separation_pixels grid_width_pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up images to be saved
pattern=abs(pattern-1);
 fig=figure;
 imshow(pattern);
 if speckles==1;
 name=['region / um: ', num2str(region_x), '     stitchHFW / um: ',num2str(hfw), '     subset size / um: ',num2str(subset_dimension)];
 else
 name=['region / um: ', num2str(region_x)];
 end
 title(name);
 
% name of generated file
hfw_round=round(hfw);
if speckles==1;
name=['reg',num2str(region_x),'_hfw',num2str(hfw_round),'_sw',num2str(speckle_width),'_sd',num2str(subset_dimension)];
if grids==1;
    name=strcat(name,['_g',num2str(grid_spacing)]);
end
end

if speckles==0;
    if grids==1;
    name=['reg',num2str(region_x),'_g',num2str(grid_spacing)];
    end
end
%
% Generate GDS file
gds_generator;

%3500 200nm square speckles takes about an hour...
if speckles==1;
rate=0.04*3500;%in um2 per hour
est_time=(length(final_speckle_i)*speckle_width^2)./rate
end

% save image and workspace file
imwrite(pattern,[name,'.bmp'] );
save([name,'.mat']);
movefile(['Cells\',name,'.gds'],[name,'.gds']);
movefile reg* Patterns_120319

% remove crud that also gets generated for some reason...
rmdir Cells;
delete(['Cells\',name])
delete(['Cells\',name,'_gds'])


