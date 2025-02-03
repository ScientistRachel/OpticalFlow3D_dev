function [vx,vy,rel] = calc_flow2D(images ,xySig, tSig, wSig)
% The calc_flow3D function calculates optical flow velocities for a single
% z-stack of images. Surrounding z-stacks in time are necessary to perform
% the calculations. To peform calculations on an entire timelapse, see
% parse_flow.m
%
% USAGE: [vx,vy,vz,rel] = calc_flow3D(images ,xyzSig, tSig, wSig)
%
% INPUTS:
% images = 4D matrix with dimensions N_X, N_Y, N_Z, N_T
%           N_T should be odd as only the central timepoint will be analyzed.
%           N_T must be greater than or equal to 3*tSig+1.
% xySig  = sigma value for smoothing in all spatial dimensions. Default 3.
%           Larger values remove noise but remove spatial detail.
% tSig   = sigma value for smoothing in the temporal dimension. Default 1.
%           Larger values remove noise but remove temporal detail.
% wSig   = sigma value for Lucas-Kanade neighborhood. Default is 4.
%           Larger values include a larger neighboorhood in the
%           Lucas-Kanade constraint and will smooth over small features.
%
% OUTPUTS:
% vx     = Velocity in the x direction, reported as pixels/frame
% vy     = Velocity in the y direction, reported as pixels/frame
% rel    = Reliability, the smallest eigenvalue of (A'wA)
%           This is a measure of confidence in the linear algebra solution
%           and can be used to mask the velocities for downstream analysis.
%
% Change Log:
% 2025/02/03 Rachel M. Lee - Created function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the images are 2D + time
if length(size(images)) ~= 3
    error('Error: Input image must be a 3D matrix with dimensions N_X, N_Y, N_T')
end

% Implement default values if they are not specified
if ~exist('xySig','var') || isempty(xySig)
    disp('Using default xyzSig value of 3...')
    xySig = 3;
end
if ~exist('tSig','var') || isempty(tSig)
    tSig = 1;
    disp('Using default tSig value of 1...')
end
if ~exist('wSig','var') || isempty(wSig)
    wSig = 4;
    disp('Using default wSig value of 4...')
end

% Check image size against tSig
Nt =  size(images,3);
if Nt < 3*tSig+1
    error('Error: Input images will lead to edge effects. N_T must be >= 3*tSig+1')
end
if ~mod(Nt,2) % enforce odd number
    error('Error: Input images must have an odd number of timepoints. Only the central time point is analyzed')
end
NtSlice = ceil(Nt/2);


%% Set up filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Common spatial terms
x = -ceil(3*xySig):ceil(3*xySig);
xySig2 = xySig/4;
y = -ceil(3*xySig2):ceil(3*xySig2);
fderiv = exp(-x.*x/2/xySig/xySig)/sqrt(2*pi)/xySig;
fsmooth = exp(-y.*y/2/xySig2/xySig2)/sqrt(2*pi)/xySig2;
gderiv = x/xySig/xySig;
gsmooth = 1;

%   Build x-gradient filter kernels (along first dimension)
xFil1 = (fderiv.*gderiv)';
yFil1 = (fsmooth.*gsmooth);
%   Build y-gradient filter kernels (along second dimension)
xFil2 = (fsmooth.*gsmooth)';
yFil2 = (fderiv.*gderiv);

%   Build t-gradient filter kernels (t = third dimension)
t = -ceil(3*tSig):ceil(3*tSig);
fx = exp(-x.*x/2/xySig/xySig)/sqrt(2*pi)/xySig;
ft = exp(-t.*t/2/tSig/tSig)/sqrt(2*pi)/tSig;
gx = 1;
gt = t/tSig/tSig;
xFil3 = (fx.*gx)';
yFil3 = xFil3';
tFil3 = permute(ft.*gt,[3 1 2]);

% Structure tensor -- Lucas Kanade neighborhood filter
wRange = -ceil(3*wSig):ceil(3*wSig);
gw = exp(-wRange.*wRange/2/wSig/wSig)/sqrt(2*pi)/wSig;
xFil4 = (gw)';
yFil4 = (gw);

clear gderiv gsmooth gt gw gx ft fx fsmooth fderiv x y t gw wRange

%% Spatial and Temporal Gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial gradients require at least N_T = 3 to avoid edge effets, while 
% the temporal gradient requires N_T >= 3*tSig+1.
% After gradients are calculated, only the central 3 slices are kept

% dx
dxI = imfilter(imfilter(images, xFil1, 'replicate'), yFil1, 'replicate'); % Filtering to calculate the gradient
dxI = dxI(:,:,NtSlice-1:NtSlice+1); % Only need the slice of interest plus one timepoint before and after for remaining calculations, simplify to save memory
clear xFil1 yFil1 zFil1

% dy
dyI = imfilter(imfilter(images, xFil2, 'replicate'), yFil2, 'replicate');
dyI = dyI(:,:,NtSlice-1:NtSlice+1);
clear xFil2 yFil2 zFil2

% dz
dtI = imfilter(imfilter(imfilter(images, xFil3, 'replicate'), yFil3, 'replicate'), tFil3, 'replicate');
clear images % No longer needed for calculations; clear to save memory
dtI = dtI(:,:,NtSlice-1:NtSlice+1);
clear xFil3 yFil3 zFil4 tFil3

%% Structure Tensor Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following calculations are for the individual elements of the
% matrices required for the optical flow calculation, incorporating
% Gaussian weighting into the Lucas-Kanade constraint.
%
% Once filtering is done, only the central timepoint is kept to save save
% memory. Because all outputs of the above section are N_X x N_Y x N_Z x 3,
% the central slice = 2.

% Time componenents
wdtx = imfilter(imfilter(dxI.*dtI, xFil4, 'replicate'), yFil4, 'replicate');
wdtx = wdtx(:,:,2);

wdty = imfilter(imfilter(dyI.*dtI, xFil4, 'replicate'), yFil4, 'replicate');
wdty = wdty(:,:,2);

clear dtI

% Spatial Components
wdxy = imfilter(imfilter(dxI.*dyI, xFil4, 'replicate'), yFil4, 'replicate');
wdxy = wdxy(:,:,2);

wdx2 = imfilter(imfilter(dxI.*dxI, xFil4, 'replicate'), yFil4, 'replicate');
wdx2 = wdx2(:,:,2);

clear dxI

wdy2 = imfilter(imfilter(dyI.*dyI, xFil4, 'replicate'), yFil4, 'replicate');
wdy2 = wdy2(:,:,2);

clear dyI
clear xFil4 yFil4


%% Optical Flow calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation is v = (A' w A)^-1 A' w b
% A = -[dxI dyI]
% b = [dtI]
% w multiplication is incorporated in the structure tensor inputs above
% A' w b = -[wdtx wdty]  (minus sign because of negative sign on A)
% (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]

% Determinant
determinant = (wdx2.*wdy2) - (wdxy.^2);

% A^-1 = [a b ; c d]^-1 = (1/det(A))[d - b; -c a]
vx = ((determinant + eps).^-1).*((wdy2.*-wdtx)+(-wdxy.*-wdty));
vy = ((determinant + eps).^-1).*((-wdxy.*-wdtx)+(wdx2.*-wdty));

clear wdtx wdty

%% Eigenvalues for Reliability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve det(A^T w A - lamda I) = 0
% (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]

trace = wdx2 + wdy2;

clear wdx2 wdxy wdxy wdy2

L1 = (trace + sqrt(trace.^2 - 4*determinant))/2;
L2 = (trace - sqrt(trace.^2 - 4*determinant))/2;
rel = real(min(L1, L2));

clear L1 L2
