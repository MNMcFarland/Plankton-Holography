function [S,D] = recsplatgpu3(fin,nr,res,wvl,w,A)

% Numerical reconstruction of in-line holograms
% Reconstructs images at distances indicated by nr
% Computes the extended depth of field image from reconstructed planes
% using the local coefficient of variation of the gradient as a focus metric
%
% USAGE: [S,D] = recsplatgpu3(fin,nr,res,wvl,w,A)
%
% INPUT:
%   fin - hologram image or file name(s) (string, cell, or dir structure)
%   nr - row vector of actual (in air) plane depths in um
%   res - resolution in um/pixel, scalar or 2 element vector
%       of resolution along y (rows) and x (columns)
%   wvl - laser wavelength in um
%   w - (optional) plane processing window radius, default = 32
%   A - (optional) image for background subtraction
%
% OUTPUT: writes reconstructed splat image to subfolder rec3 if no output arguments
%   S - extended depth of focus splat image
%   D - depth map, plane depth for each pixel

if ~exist('w','var') || isempty(w); w = 32; end % plane processing window radius
if ~exist('A','var') || isempty(A); A = 0; end % 
A = A - mean(A,'all'); % zero center background image
A = gpuArray(A);

if isnumeric(fin) && all(size(fin)>1) 
    I = fin;
    fnm = 0;
    imx = max(I(:));
elseif isstruct(fin) && isfield(fin,'name')
    fnm = {fin.name};
elseif ~iscell(fin)
    fnm = {fin};
else
    fnm = fin;
end

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
    mkdir('rec3'); % create rec directory if it doesn't exist
end

% set up a progress bar with cancel button
wb = waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','generating splat image');
hnd = findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
cntr = 0; % counter for progress bar
tot = length(r)*length(fnm); % total # of planes to reconstruct

% read first image
if ~exist('I','var') 
    I = imread(fnm{1});
    imx = double(intmax(class(I)));
end

% generate the reconstruction filter kernel
Nx = size(I,2); % Image Size # columns
Ny = size(I,1); % Image Size # rows
x = ((1:Nx)-Nx/2)/(Nx*dx);
y = ((1:Ny)-Ny/2)/(Ny*dy);
[x,y] = meshgrid(x,y);
arg = -1i*wvl*pi*(x.^2+y.^2); % for Kirchhoff-Fresnel kernel
arg = fftshift(arg);
arg = gpuArray(single(arg));

% for Rayleigh-Sommerfeld kernel only
% k=gpuArray(2*pi/wvl);
% arg=sqrt(1-(wvl*x).^2-(wvl*y).^2);

nbhd = gpuArray(ones(w*2+1));
mnflt = gpuArray(nbhd/(w*2+1)^2); % kernel to calculate local mean with imfilter

for f=1:length(fnm)
    if f>1; I = imread(fnm{f}); end % first image already open
    I = single(I);
    I = gpuArray(I);
    I = I - A; % background subtraction

    D = zeros(Ny,Nx); % initialize depth map
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
%         M1 = stdfilt(M1,nbhd) ./ imfilter(M1,mnflt,'symmetric'); % CV of gradient (Tamura coefficient of the gradient (squared))
        % stdev(X) = sqrt(E[X^2] - E[X]^2); % alternate formulation of the stdev
%         M1 = M1 - 1; % subtract constant k = 1 for numerical stability, to avoid catastrophic cancellation
        Eg = imfilter(G,mnflt,'symmetric'); % local mean of gradient, E[G]
        Eg2 = imfilter(G.^2,mnflt,'symmetric'); % local mean of squared gradient, E[G^2]
        % G = sqrt(Eg2 - Eg.^2)./Eg; % local coefficient of variation of gradient = std/mean
        G = (Eg2 - Eg.^2)./Eg; % local variance to mean ratio of gradient = var/mean

        ind = G>F; % find max focus metric
        F(ind) = G(ind); % store max focus metric
        
        S(ind) = P(ind); % update Splat image
        D(ind) = z*n; % update the depth map
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(f) ' of ' num2str(length(fnm))...
            '; z = ' num2str(z*n)]) 
    end
    S = gather(S);
    D = gather(D);
    if nargout==0
%         S = S/4095*65535; % scale splat from 12 bit to 16 bit
        S = S/imx*65535; % scale splat from original bit depth to 16 bit
        [~,fstr] = fileparts(fnm{f});
        fstr = [fstr '_' num2str(min(nr)/1000) '-' num2str(max(nr)/1000)]; % append depth range to file name in mm
%         imwrite(uint16(S),['rec2\splat2_' fstr '.tif'],'tif','compression','none')
%         imwrite(uint16(D),['rec2\splat2D_' fstr '.tif'],'tif','compression','none')
        imwrite(uint16(S),['rec3\splat3_' fstr '.png'],'png')
        imwrite(uint16(D),['rec3\splat3D_' fstr '.png'],'png') % limited to 6.55 cm max pathlength unless scaled first
    end
end
delete(wb)
warning('on','MATLAB:MKDIR:DirectoryExists')
