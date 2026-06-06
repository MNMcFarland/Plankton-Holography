function [S1i,S2i,D]=recsplat(fin,nr,res,lambda,m,t,w)

% Numerical reconstruction of holograms
% Reconstructs images at depths indicated by nr
% Computes the splat image of reconstructed planes
%
% USAGE: [S1,S2,,D]=recsplat(fin,nr,res,lambda,m,t,w)
%
% INPUT:
%   fin - hologram image file name(s) (string, cell, or dir structure)
%   nr - row vector of plane depths in um
%   res - resolution in um/pixel, scalar or 2 element vector
%       of resolution along y (rows) and x (columns)
%   lambda - laser wavelength in um
%   m - (optional) plane processing flag,  default = true
%   t - (optional) plane processing threshold, default = 0.091554*intmax(fin)
%   w - (optional) plane processing window size, default = 32
%
% OUTPUT: writes reconstructed splat image to subfolder rec
%   S1 - minimum pixels splat image
%   S2 - max stdev above threshold splat image
%   D - depth map

if ~exist('m','var') || isempty(m); m=1; end
if (~exist('w','var') || isempty(w)) && m; w=32; end % plane processing window size

if isstruct(fin) && isfield(fin,'name')
    fnm={fin.name};
elseif ~iscell(fin)
    fnm={fin};
else
    fnm=fin;
end

n=1.33; % refractive index of medium (i.e. water)
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
% mkdir('rec_w'); % create rec_w directory if it doesn't exist

% set up a progress bar with cancel button
wb=waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','reconstructing hologram');
hnd=findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
cntr=0; % counter for progress bar
tot=length(r)*length(fnm); % total # of planes to reconstruct

I=imread(fnm{1});
[Ny,Nx]=size(I); % Ny = # rows; Nx = # columns
x =((1:Nx)-Nx/2)/(Nx*dx);
y =((1:Ny)-Ny/2)/(Ny*dy);
[x,y]=meshgrid(x,y);
arg=-1i*lambda*pi*(x.^2+y.^2); % for Kichhoff-Fresnel kernel

% for Rayleigh-Sommerfeld kernel only
% k=2*pi/lambda;
% arg=sqrt(1-(lambda*x).^2-(lambda*y).^2);

imx=double(intmax(class(I)));
if (~exist('t','var') || isempty(t)) && m; t=.1221*imx; end % threshold for reconstructed plane processing = 8000 for uint16
% if (~exist('t','var') || isempty(t)) && m; t=.091554*imx; end % threshold for reconstructed plane processing = 6000 for uint16

% bpfun = @(bs) ones(size(bs.data)).*sum(bs.data(:)>t); % block processing function
% M2 = zeros(Ny,Nx);

for f=1:length(fnm)
    if f>1; I=imread(fnm{f}); end
    [~,fstr]=fileparts(fnm{f});
    fstr=[fstr '_' num2str(min(nr)/1000) '-' num2str(max(nr)/1000)]; % depth range in mm
    
%     I = double(I);
    I = fft2(I);
    
    S1 = ones(Ny,Nx)*imx; % initialize splat image
    D = zeros(Ny,Nx); % initialize depth map
    
    if m==1
        S2 = zeros(Ny,Nx); % initialize splat image
        G = zeros(Ny,Nx); % initialize array to store pixel counts
    end
    
    for z = r
        if getappdata(wb,'canceling') % check if cancelled
            S1i=0.75*imx-S1;
            if m==1
                S2i = .75*imx-S2; % invert splat image
            end
            delete(wb)
            error('cancelled')
%             return
        end

%         n = (2*pi*1i/k)*exp(-1i*2*pi^2*z*(x.^2 + y.^2)/k); % original
%         n = 1i*lambda*exp(-1i*lambda*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel for phase
%         n = exp(-1i*lambda*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel
%         n = exp(-1i*k*z*arg); % Rayleigh-Sommerfeld kernel (needs arg defined above)
%         n = ifftshift(n);
%         M = abs(ifft2(I.*n));
%         M = abs(atan(imag(M)./real(M))); % phase (only seems to work with K-F kernel)

        K = exp(z*arg); % Kirchhoff-Fresnel kernel
%         K=exp(-1i*k*z*arg); % Rayleigh-Sommerfeld kernel
        K = ifftshift(K);
        K = single(K);
        M = I.*K;
        M = ifft2(M);
%         M=abs(atan(imag(M)./real(M))); % phase? (only seems to work with K-F kernel)
        M = abs(M);

        %% update splat image
        ind1 = M<S1; % find min
        S1(ind1) = M(ind1); % update splat image

        %% process reconstructed image
        if m==1
            M1 = imgaussfilt(M,4); % smooth
            M1 = imgradient(M1); % gradient magnitude
            M1 = imfilter(single(M1>t), ones(w), 'symmetric'); % count pixels > threshold within local nhood (w x w)
%             M1 = blockproc(M1, [w w], bpfun); % count pixels > threshold within local nhood
%             for u = 1:w:Ny-w
%                 for v = 1:w:Nx-w
%                     M2(u:u+w-1,v:v+w-1) = sum(M1(u:u+w-1,v:v+w-1)>t,'all');
%                 end
%             end
%             M1 = M2;

            ind2 = M1>G;
            G(ind2) = M1(ind2); % update obj func max
            
            S2(ind2) = M(ind2); % update Splat image
            D(ind2) = z*n; % update depth map
        else
            D(ind1) = z*n;
        end

        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(f) ' of ' num2str(length(fnm))]) 
    end
    S1i = .75*imx-S1; % invert and scale splat image
    S1i = uint16(S1i);
    if m==1
        S2i = .75*imx-S2; % invert splat image
        S2i = uint16(S2i);
    else
        S2i = [];
    end
    if nargout==0
%         imwrite(S1i,['rec\splat1_' fstr '.jpg'],'jpg','quality',90)
        imwrite(S1i,['rec\splat1_' fstr '.tif'],'tif','compression','none')
        imwrite(uint16(D),['rec\splatD_' fstr '.tif'],'tif','compression','none')
        if exist('S2i','var') && ~isempty(S2i)
            imwrite(S2i,['rec\splat2_' fstr '.tif'],'tif','compression','none')
        end
    end
end
delete(wb)
warning('on','MATLAB:MKDIR:DirectoryExists')
