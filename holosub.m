function A = holosub(fnm,bkg,num)

% compute and subtract time averaged background
% when called with fnm argument: writes background subtracted images to \sub directory
% when called with empty fnm argument: returns average image
%
% USAGE: A = holosub1(fnm,bkg,num)
%
% EXAMPLES:
%   A = holosub1([],bkg,num); compute background from files in bkg and return average A
%   holosub1(fnm); compute background, subtract from each file in fnm, and save to /sub1 directory
%   holosub1(fnm,bkg,num); specify files to use to generate background
%
% INPUT:
%   fnm - (optional) name(s) of files to background subtract (cell or dir structure)
%   bkg - (optional) names of files used to compute background (cell or dir structure)
%   num - (optional) maximum number of frames to average (default=50 if bkg not specified)
%
% OUTPUT:
%   with fnm input argument: background subtracted images written to subfolder 'sub'
%   otherwise:
%   A - average image

if isstruct(fnm) && isfield(fnm,'name')
    fnm = {fnm.name};
elseif ~iscell(fnm) && ~isempty(fnm)
    fnm = {fnm}; 
end

if ~exist('bkg','var') || isempty(bkg)
    bkg = fnm;
    num = 50;
elseif isstruct(bkg) && isfield(bkg,'name')
    bkg = {bkg.name};
elseif ~iscell(bkg)
    bkg = {bkg}; 
end

if ~exist('num','var') || num>length(bkg)
    num = length(bkg); 
end

% set up a progress bar with cancel button
wb = waitbar(0,'computing background','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','processing holograms');
hnd = findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
cntr = 0; % counter for progress bar
tot = length(fnm)+num; % total # of images to process

%% compute average background

A = double(imread(bkg{1}));
for n = round(linspace(2,length(bkg),num-1))
    if getappdata(wb,'canceling') % check if cancelled
        delete(wb)
        return
    end
    I = double(imread(bkg{n}));
    A = A + I;
    
    % update progress bar
    cntr = cntr+1;
    waitbar(cntr/tot,wb,['computing background']) 
end

A = A/num; % average image
% G = mean(A,'all');

%% 3. subtract average and write image files to sub directory
if ~isempty(fnm)
    Az = A - mean(A,'all'); % zero centered average image
    mkdir('sub');
    for n = 1:length(fnm)
        if getappdata(wb,'canceling') % check if cancelled
            delete(wb)
            return
        end

        I = double(imread(fnm{n}));
        I = I - Az;
        I = uint8(I);
%         imwrite(I,['sub1\' fnm{n}(1:end-4) '.tif'],'tif','Compression','none');
        imwrite(I,['sub' filesep fnm{n}(1:end-4) '.bmp'],'bmp');
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,'subtracting background') 
    end
end
delete(wb)
