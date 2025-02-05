function ExportMap()
%%EXPORTMAP converts the master GDS into a specific fab GDS
% Author : Nicolas Ayotte                                  Creation date : 01/04/2014
% 
% The entire programs is in microns.


%% Initialize Cell Parameters
%tp = cd; cd ../Functions/; addpath(genpath(cd)); cd(tp);     % add path to every folder of the library

tic;
log = SetupLog('do', true);  % options 'do' and 'file'; functions 'write' and 'close'
log.write('\n%s\nCELL FUNCTION %s\n\n', log.bar(), log.title());
cad = ProjectDefinition(log);


%% Merge all the GDS cells
InitializeGDSclasses();
CastLayerMap(cad.outfil, cad.layermap, 'output', log);


%% Close the log file
log.write('\nEND  -  %s\n\n', log.time());
log.close();

end
