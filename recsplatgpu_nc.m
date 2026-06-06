function [S,D] = recsplatgpu_nc(fin,nr,res,wvl,w)

% Numerical reconstruction of in-line holograms stored in .nc file
% Reconstructs images for each frame at depths indicated by nr
% Computes the splat image of reconstructed planes for all frames
% using the local coefficient of variation of the gradient as a focus metric
%
% USAGE: [S,D] = recsplatgpu_nc(fin,nr,res,wvl,w)
%
% INPUT:
%   fin - single netcdf format (.nc) file name (char, string, or dir structure)
%   nr - row vector of actual (in air) plane depths in um
%   res - resolution in um/pixel, scalar or 2 element vector
%       of resolution along y (rows) and x (columns)
%   wvl - laser wavelength in um
%   w - (optional) plane processing window radius, default = 32
%
% OUTPUT: writes reconstructed splat image to subfolder rec
%   S - extended depth of focus splat image
%   D - depth map, plane depth for each pixel

if ~exist('w','var') || isempty(w); w = 32; end % plane processing window radius

if isstruct(fin) && isfield(fin,'name') && isfield(fin,'folder')
    fnm = [fin.folder filesep fin.name];
elseif ischar(fin) %~iscell(fin)
    fnm = fin;
elseif isstring(fin) && length(fin)==1
    fnm = char(fin);
else
    error('recsplatgpu_nc : unrecognized input');
end
nci = ncinfo(fnm);
dim = [nci.Groups.Dimensions.Length]; % [frames width height ctd]
nfrm = dim(1);
[~,dnm] = fileparts(fnm); % directory name

% if ~exist('A','var') || isempty(A); A = zeros(dim(3),dim(2)); end
A = holosub_nc([],fnm);
A = single(A-mean(A,'all')); % zero centered mean

% refractive index of medium alters apparent plane depth
n = 1.333; % refractive index of medium (water)
r = nr/n; % divide by refractive index of medium

if length(res)==2
    dy = res(1);  % Pixel resolution along Y (rows)
    dx = res(2);  % Pixel resolution along X (columns)
elseif length(res)==1
    dy = res;  % Pixel resolution along Y
    dx = res;  % Pixel resolution along X
else
    error('res must be scalar or a 2 element vector')
end

if nargout==0
    warning('off','MATLAB:MKDIR:DirectoryExists')
    mkdir(dnm); % create rec directory if it doesn't exist
end

% set up a progress bar with cancel button
wb = waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','generating splat image');
hnd = findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
cntr = 0; % counter for progress bar
tot = length(r)*nfrm; % total # of planes to reconstruct

% generate the reconstruction filter kernel
Nx = dim(2); % size(I,2); Image Size # columns
Ny = dim(3); % size(I,1); Image Size # rows
x = ((1:Nx)-Nx/2)/(Nx*dx);
y = ((1:Ny)-Ny/2)/(Ny*dy);
[x,y] = meshgrid(x,y);
arg = -1i*wvl*pi*(x.^2+y.^2); % for Kirchhoff-Fresnel kernel
arg = fftshift(arg);
arg = gpuArray(single(arg));

% for Rayleigh-Sommerfeld kernel only
% k=gpuArray(2*pi/wvl);
% arg=sqrt(1-(wvl*x).^2-(wvl*y).^2);

nbhd = gpuArray(true(w*2+1));
mnflt = gpuArray(nbhd/(w*2+1)^2); % kernel to calculate local mean with imfilter

for f = 1:nfrm
    I = ncread(fnm,[nci.Groups.Name '/Images'],[1 1 f],[dim(2) dim(3) 1])'; 
    I = single(I) - A; % subtract background
    I = gpuArray(I);

    D = gpuArray.zeros(Ny,Nx); % initialize depth map
%     S=gpuArray.zeros(Ny,Nx); % initialize splat image
    S = I; % initialize splat image
    F = gpuArray.zeros(Ny,Nx,'single'); % initialize array to store focus metric
    
    I = fft2(I);
    
    for z = r
        if getappdata(wb,'canceling') % check if cancelled
            S = gather(S);
            D = gather(D);
            delete(wb)
%             error('cancelled')
            return
        end
        
%         K = (2*pi*1i/k)*exp(-1i*2*pi^2*z*(x.^2 + y.^2)/k); % original
%         K = wvl*1i*exp(-1i*lambda*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel for phase
%         K = exp(-1i*wvl*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel
        K = exp(z*arg); % Kirchhoff-Fresnel kernel
%         K = exp(-1i*k*z*arg); % Rayleigh-Sommerfeld kernel
%         K = fftshift(K);
        K = single(K);
        P = I.*K;
        
        P = ifft2(P);
        P = real(P);

        % compute focus metric
        G = imgradient(P,'central'); % gradient magnitude
%         G = stdfilt(G,nbhd) ./ imfilter(G,mnflt,'symmetric'); % CV of gradient (Tamura coefficient of the gradient (squared))
        Eg = imfilter(G,mnflt,'symmetric'); % local mean of gradient, E[G]
        Eg2 = imfilter(G.^2,mnflt,'symmetric'); % local mean of squared gradient, E[G^2]
        G = sqrt(Eg2 - Eg.^2)./Eg; % local coefficient of variation of gradient = std/mean

        ind = G>F; % find max focus metric
        F(ind) = G(ind); % store max focus metric
        
        S(ind) = P(ind); % update Splat image
        D(ind) = z*n; % update the depth map
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(f) ' of ' num2str(nfrm)...
            '; z = ' num2str(z*n)]) 
    end
    S = gather(S);
    D = gather(D);
    if nargout==0
        S = S/4095*65535; % scale splat from 12 bit to 16 bit
        [~,fstr] = fileparts(fnm);
        fstr = [fstr '_' num2str(f) '_' num2str(min(nr)/1000) '-' num2str(max(nr)/1000)]; % append depth range to file name in mm
        imwrite(uint16(S),[dnm '\splat3_' fstr '.png'],'png')
        imwrite(uint16(D),[dnm '\splat3D_' fstr '.png'],'png') % limited to 6.55 cm max pathlength unless scaled first
    end
end
delete(wb)
warning('on','MATLAB:MKDIR:DirectoryExists')
