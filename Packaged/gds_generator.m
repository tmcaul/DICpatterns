%% Generate GDS file format of grids and speckles
% Author : Tom McAuliffe
% Creation date : 06/06/2018
% The entire programs is in microns.

%function takes in dimensions of patterning region, speckle locations,
%whether grids are required, and if so the spacing and width of the grids (in um).

%function completed=gds_generator(region_x, region_y, speckle_x, speckle_y, grids, grid_spacing, grid_width);

% add path to every folder of the library
[log, cad, cellname, topcell, layerMap] = InitializeCell(); %draws on project definition, which should be in the same folder as this.

% Load GDS Library References
refs = [];
refs = GetRefsFloorplan(refs);

% Start making layers

f=waitbar(1./(length(speckle_x)+2.*(round(region_x./grid_spacing))), 'generating pattern...');
a=0;
b=0;
c=0;

if speckles==1;
%generates a speckle pattern from final_speckle_x, final_speckle_y
for a=1:length(speckle_x);
waitbar(a./(length(speckle_x)+2.*(round(region_x./grid_spacing))));
    
infoIn=CursorInfo([speckle_x(a),speckle_y(a)],0,1);
infoIn.ori=-45;
[topcell,infoOut,infoIn]=PlaceRect(topcell,infoIn,speckle_width,speckle_width,0,0);
a=a+1;
end
end

%generate grids
if grids==1;
for b=1:floor(region_x./grid_spacing);
waitbar((a+b)./(length(speckle_x)+2.*(round(region_x./grid_spacing))));
infoIn=CursorInfo([b.*grid_spacing,region_y./2],0,4);
[topcell,infoOut,infoIn]=PlaceRect(topcell,infoIn,grid_width,region_y,0,0);
b=b+1;
end

for c=1:floor(region_y./grid_spacing);
waitbar((a+b+c)./(length(speckle_x)+2.*(round(region_x./grid_spacing))));
infoIn=CursorInfo([0,c.*grid_spacing],0,1);
[topcell,infoOut,infoIn]=PlaceRect(topcell,infoIn,region_x,grid_width,0,0);
c=c+1;
end
end
close(f)

% Save GDS and .mat cell information
cellname=name;
topcell.sname=name;
FinalizeCell(cad, cellname, topcell, refs, infoIn, infoOut, log);
%MergeCells();
%ExportMap();
completed=1;
%end