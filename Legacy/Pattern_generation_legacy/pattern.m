clear all

% horizontal field width in um
hfw=50;
img_size=[884, 1024]; %rows and columns so x, y reversed

%want 3x3 speckles per subset, 3x3 pixels per subset
%square subset dimension in um
subset_dimension=0.5;

no=2; %no of speckles (-1) that affect shifts (2 works well)
thresh=0.5; %of +/- of equilibrium separation that will shift 
mag=0.15; %of shift of micro-adjustments


%
d_um=9./subset_dimension; %density of speckles: no. per um^2
pixel_width=(1/img_size(2))*hfw; % width of one pixel in um (um per pixel)
speckle_width=3.*pixel_width; % in um

%subset size determined by speckle spacing.

d_pixels=d_um.*(pixel_width.^2); %density of speckles in no. per sq. pixel.
speckle_separation=round(1/d_pixels); %speckle separation in pixels

%each speckle should cover at least 3x3 pixels
%each subset should contain at least 3x3 speckles

%
% set up lattice of speckle centres with density set by d (and therefore speckle separation)
% each speckle then can move by a random number of pixels in x and y
% amount by which each pixel can move is 1/3 of ordered distance

image=zeros(img_size);
c=1; %initialise offset trackers

for i=1:img_size(1);
    for j=1:img_size(2);
        
        % random offset between +/- 1/3 of speckle separation
        if mod(i,speckle_separation)==0 & mod(j,speckle_separation)==0;
            
            offset1(c)=randi([-round(speckle_separation/3), round(speckle_separation/3)]);
            offset2(c)=randi([-round(speckle_separation/3), round(speckle_separation/3)]);
            
            image(i+offset1(c)-round(speckle_separation/3),j+offset2(c)-round(speckle_separation/3))=1; 
        end
        j=j+1;
    end
    i=i+1;
end
clear offset1 offset2 i j a speckles_init speckles
%

%
original_image=image;

[speckle_i,speckle_j]=find(image); %find non-zero elements and their indices
speckle_locs(:,1)=speckle_i;
speckle_locs(:,2)=speckle_j;
clear speckle_i speckle_j

% initialise variables (so code repeatable)
distances=[0,0];
distances_sorted_cropped=[0,0];
speckles_closest=[0,0];
sq_distances=[0,0];
vectors=[0,0];
net_vector=[0,0];
loc=0;
speckle_pos_new=[0,0];

% randomise speckle to optimise

%%%% do micro-adjustments on random speckles
for index=1:30000;
    
    a=randi([1,length(speckle_locs)]);
    
    loc=speckle_locs(a,:); %considered speckle
    distances=speckle_locs-loc; %distances
    
    sq_distances=sqrt((distances(:,1)).^2+(distances(:,2)).^2); %square distances
    [distances_sorted(:,1),distances_sorted(:,2)]=sort(sq_distances); %outputs min. square distances and the index in the array
    
    distances_sorted_cropped=distances_sorted(1:no,:); %crops array to the speckle itself and no. closest
    
    for b=1:no; %closest speckles i,j indices
        speckles_closest(b,:)=speckle_locs(distances_sorted(b,2),:);
    end
    
    vectors=speckles_closest-speckles_closest(1,:); %get position vectors, first row is always speckle itself
    net_vector=sum(vectors(2:end,:));
    
    if norm(net_vector)>(1+thresh)*speckle_separation; %if net vector magnitude > (1+thresh) * desired separation  
        speckle_pos_new=round(loc+mag*net_vector); % slight shift in direction of net vector (to move it closer)
        image(speckle_pos_new(1),speckle_pos_new(2))=1;% create speckle at new location
        image(loc(1),loc(2))=0; %remove old speckle    
    end
    
    if norm(net_vector)<(1-thresh)*speckle_separation; %if other speckles are too close...
        speckle_pos_new=round(loc-mag*net_vector);
        image(speckle_pos_new(1),speckle_pos_new(2))=1;% create speckle at new location
        image(loc(1),loc(2))=0; %remove old speckle    
    end
    
end

clear a b speckles_closest sq_distances loc speckle_locs speckle_pos_new c distances distances_sorted distances_sorted_cropped net_vector vectors
%
% convolution of 3x3 pixels speckle at each centre (which is optimal)
pixel=ones(3,3);
pattern=conv2(pixel,image);
pattern_1=conv2(pixel,original_image);

pattern=pattern(1:img_size(1),1:img_size(2));
pattern_1=pattern_1(1:img_size(1),1:img_size(2));

clear pixel

%
figure;
%subplot(1,3,1)
imshow(pattern);
title(['size: ', num2str(img_size), '     HFW / um: ',num2str(hfw),'     speckle width / um: ', num2str(speckle_width), '     pixel width / um: ',num2str(pixel_width)]);
%subplot(1,3,2)
figure;
imshow(pattern_1);
title('initial');
%subplot(1,3,3)
figure;
imshow(pattern-pattern_1);
title('difference');
