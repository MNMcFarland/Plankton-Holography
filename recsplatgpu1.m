function [Si,D]=recsplatgpu1(fin,nr,res,wvl)

% Numerical reconstruction of holograms
% Reconstructs images at depths indicated by nr
% Computes the splat image of reconstructed planes by
% finding the minimum intensity for each pixel in all planes.
%
% USAGE: [S1,S2,D]=recsplatgpu(fin,nr,res,lambda,m,t,w)
%
% INPUT:
%   fin - hologram image file name(s) (string, cell, or dir structure)
%   nr - row vector of actual (in air) plane depths in um
%   res - resolution in um/pixel, scalar or 2 element vector
%       of resolution along y (rows) and x (columns)
%   lambda - laser wavelength in um
%   m - (optional) plane processing flag,  default = true
%   t - (optional) plane processing threshold, default = 0.091554*intmax(fin)
%   w - (optional) plane processing window size, default = 32
%
% OUTPUT: writes reconstructed splat image to subfolder rec
%   S1 - per pixel minimum splat image
%   S2 - max stdev above threshold splat image
%   D - depth map, plane depth for each pixel

if isstruct(fin) && isfield(fin,'name')
    fnm={fin.name};
elseif ~iscell(fin)
    fnm={fin};
else
    fnm=fin;
end

% refractive index of medium alters apparent plane depth
n=1.333; % refractive index of medium (water)
r=nr/n; % divide by refractive index of medium

if length(res)==2
    dy = res(1);  % Pixel resolution along Y (rows)
    dx = res(2);  % Pixel resolution along X (columns)
elseif length(res)==1
    dy = res;  % Pixel resolution along Y
    dx = res;  % Pixel resolution along X
else
    error('res must be scalar or a 2 element vector')
end

warning('off','MATLAB:MKDIR:DirectoryExists')
mkdir('rec'); % create rec directory if it doesn't exist

% set up a progress bar with cancel button
wb=waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','reconstructing hologram');
hnd=findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
cntr=0; % counter for progress bar
tot=length(r)*length(fnm); % total # of planes to reconstruct

I=imread(fnm{1});
Nx=size(I,2); % Image Size # columns
Ny=size(I,1); % Image Size # rows
x=gpuArray(((1:Nx)-Nx/2)/(Nx*dx));
y=gpuArray(((1:Ny)-Ny/2)/(Ny*dy));
[x,y]=meshgrid(x,y);
arg=-1i*wvl*pi*(x.^2+y.^2); % for Kichhoff-Fresnel kernel

% for Rayleigh-Sommerfeld kernel only
% k=gpuArray(2*pi/lambda);
% arg=sqrt(1-(lambda*x).^2-(lambda*y).^2);

imx=double(intmax(class(I)));

for f=1:length(fnm)
    if f>1; I=imread(fnm{f}); end
    [~,fstr]=fileparts(fnm{f});
    fstr=[fstr '_' num2str(max(nr)/1000-min(nr)/1000)]; % append depth range to file name
%     fstr=[fstr(1:9) fstr(end-6:end) '_' num2str(min(nr)/1000) '-' num2str(max(nr)/1000)]; % depth range
%     idx=regexp(fstr,'_IMPERX');
%     fstr=[fstr(1:idx) fstr(end-6:end) '_' num2str(max(nr)/1000-min(nr)/1000)]; % depth range

    I=single(I);
    I=gpuArray(I);
    I=fft2(I);
    
    S=gpuArray.ones(Ny,Nx)*imx; % initialize splat image
    D=zeros(Ny,Nx); % initialize depth map
    
    for z=r
        if getappdata(wb,'canceling') % check if cancelled
            Si = .75*imx-S; % invert splat image
            Si = uint16(gather(Si));
            delete(wb)
            error('cancelled')
%             return
        end
        
%         K = (2*pi*1i/k)*exp(-1i*2*pi^2*z*(x.^2 + y.^2)/k); % original
%         K = lambda*1i*exp(-1i*lambda*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel for phase
%         K = exp(-1i*lambda*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel
        K = exp(z*arg); % Kirchhoff-Fresnel kernel
%         K = exp(-1i*k*z*arg); % Rayleigh-Sommerfeld kernel
        K = ifftshift(K);
        K = single(K);
        M = I.*K;
        M = ifft2(M);
%         M = abs(atan(imag(M)./real(M))); % phase (only seems to work with K-F kernel)
        M = abs(M);
        
        ind = M<S; % find min
        S(ind) = M(ind); % update splat image
        D(ind) = z*n; % update the depth map
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(f) ' of ' num2str(length(fnm))]) 
    end
    Si = .75*imx-S; % invert splat image
    Si = uint16(gather(Si));
    if nargout==0
        imwrite(Si,['rec\splat1_' fstr '.tif'],'tif','compression','none')
        imwrite(uint16(D),['rec\splat1D_' fstr '.tif'],'tif','compression','none')
    end
end
delete(wb)
warning('on','MATLAB:MKDIR:DirectoryExists')
