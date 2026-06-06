function [S,D] = recsplatgpu2(fin,nr,res,wvl,sd,w)

% Numerical reconstruction of in-line holograms
% Reconstructs images at depths indicated by nr
% Computes the splat image of reconstructed planes
% using the local coefficient of variation of the gradient
%
% USAGE: [S,D] = recsplatgpu2(fin,nr,res,wvl,sd,w)
%
% INPUT:
%   fin - hologram image or file name(s) (string, cell, or dir structure)
%   nr - row vector of actual (in air) plane depths in um
%   res - resolution in um/pixel, scalar or 2 element vector
%       of resolution along y (rows) and x (columns)
%   wvl - laser wavelength in um
%   sd - (optional) std devs for difference of gaussians bandpass filter in spatial
%       domain (pixels), default = [0.5 2048]
%   w - (optional) plane processing window radius, default = 32
%
% OUTPUT: writes reconstructed splat image to subfolder rec
%   S - extended depth of focus splat image
%   D - depth map, plane depth for each pixel

if ~exist('sd','var') || isempty(sd); sd = [0.5 2048]; end
if ~exist('w','var') || isempty(w); w = 32; end % plane processing window radius

if isnumeric(fin) && all(size(fin)>1) 
    I = fin;
    fnm = 0;
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
    mkdir('rec2'); % create rec directory if it doesn't exist
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
if ~exist('I','var'); I = imread(fnm{1}); end

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

% difference of gaussians bandpass filter for smoothing
x = ((1:Nx)-Nx/2)/Nx;
y = ((1:Ny)-Ny/2)/Ny;
[x,y] = meshgrid(x,y);
sd1 = 1/(2*pi*sd(1)); % standard deviation in frequency domain
G1 = exp(-(x.^2+y.^2)/(2*sd1^2));
% G1(floor(Nx/2)+1,floor(Ny/2)+1) = 0; % zero center
% Gf = gpuArray(fftshift(G1));
sd2 = 1/(2*pi*sd(2)); % standard deviation in frequency domain
G2 = exp(-(x.^2+y.^2)/(2*sd2^2));
Gf = gpuArray(fftshift(G1-G2));
clear x y

% Kc = ones(w*2+1,'gpuArray'); % kernel for counting pixels above t in gradient image
% nsz = (w*2+1)^2;
nbhd = gpuArray(ones(w*2+1));
mnflt = gpuArray(nbhd/(w*2+1)^2); % kernel to calculate local mean with imfilter
% mn = mean(I(:));

for f=1:length(fnm)
    if f>1; I = imread(fnm{f}); end % first image already open
    I = single(I);
    I = gpuArray(I);

    D = zeros(Ny,Nx); % initialize depth map
%     S=gpuArray.zeros(Ny,Nx); % initialize splat image
    S = I; % initialize splat image
    G = gpuArray.zeros(Ny,Nx,'single'); % initialize array to store pixel counts
    
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
        M = I.*K;
        
        % process reconstructed image
        M1 = M.*Gf; % apply band pass filter
        M1 = real(ifft2(M1));

        M = ifft2(M);
        M = real(M);

        % compute focus metric
%         M1(M1>0) = 0; % ignore positive intensity deviations
%         M1 = stdfilt(M1,nbhd);
%         M1 = imfilter(single(M1>t), Kc, 'symmetric'); % count pixels > threshold within local nhood (w x w)
%         M1 = imfilter(single(M1.*(M1>t)), Kc, 'symmetric'); % sum values > threshold within local nhood (w x w)

        M1 = imgradient(M1,'central'); % gradient magnitude
%         M1 = stdfilt(M1,nbhd) ./ imfilter(M1,mnflt,'symmetric'); % CV of gradient (Tamura coefficient of the gradient (squared))
        % stdev(X) = sqrt(E[X^2] - E[X]^2); % alternate formulation of the stdev
%         M1 = M1 - 1; % subtract constant k = 1 for numerical stability, to avoid catastrophic cancellation
        Ex = imfilter(M1,mnflt,'symmetric'); % local mean of gradient, E[X]
        Ex2 = imfilter(M1.^2,mnflt,'symmetric'); % local mean of squared gradient, E[X^2]
        M1 = sqrt(Ex2 - Ex.^2)./Ex; % local coefficient of variation of gradient = std/mean

        ind = M1>G; % find max focus metric
        G(ind) = M1(ind); % store max focus metric
        
        S(ind) = M(ind); % update Splat image
        D(ind) = z*n; % update the depth map
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(f) ' of ' num2str(length(fnm))...
            '; z = ' num2str(z*n)]) 
    end
    S = gather(S);
    D = gather(D);
    if nargout==0
        S = S/4095*65535; % scale splat from 12 bit to 16 bit
        [~,fstr] = fileparts(fnm{f});
        fstr = [fstr '_' num2str(min(nr)/1000) '-' num2str(max(nr)/1000)]; % append depth range to file name in mm
%         imwrite(uint16(S),['rec2\splat2_' fstr '.tif'],'tif','compression','none')
%         imwrite(uint16(D),['rec2\splat2D_' fstr '.tif'],'tif','compression','none')
        imwrite(uint16(S),['rec2\splat2_' fstr '.png'],'png')
        imwrite(uint16(D),['rec2\splat2D_' fstr '.png'],'png')
    end
end
delete(wb)
warning('on','MATLAB:MKDIR:DirectoryExists')
