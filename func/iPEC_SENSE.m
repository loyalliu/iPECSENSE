function img = iPEC_SENSE(pseudoPi_kxkyc_full, adjustedCSM_xyc_full, nIter, wavWeight, varargin)

% Normalize coil sensitivity maps
tmp = sos(adjustedCSM_xyc_full);
tmp = tmp(tmp~=0);
csm = adjustedCSM_xyc_full/(mean(tmp(:)));

SensMap = csm;
RecMtx = conj(csm);

if nargin > 4
    img = varargin{1}; % Use a pre-defined initialization
else
    img = SENSE_for_POCSENSE(pseudoPi_kxkyc_full, RecMtx);
end
% find the closest diadic size for the images
[sx,sy,ncc] = size(pseudoPi_kxkyc_full);
ssx = 2^ceil(log2(sx));
ssy = 2^ceil(log2(sy));
ss = max(ssx, ssy);
W = Wavelet('Daubechies',4,4);
%W = Wavelet('Haar',2,3);
% wavWeight = 0.0001;


for iter = 1:nIter
%     disp(['iPEC-SENSE Iteration: ', num2str(iter)]);
    if wavWeight > 0
        img = zpad(img,ss,ss, 2); % zpad to the closest diadic
        img = W*(img); % apply wavelet
        img = SoftThresh(img,wavWeight); % threshold ( joint sparsity)
        img = W'*(img); % get back the image
        img = crop(img, sx, sy, 1); % return to the original size
    end
    
    % Data regeneration
    img = repmat(img, [1 1 ncc]).*SensMap;
    ksp= fft2c(img);
    
    % Data consistency
    mask = abs(pseudoPi_kxkyc_full) > 0;
    ksp(mask) = pseudoPi_kxkyc_full(mask);
    
    % Coil combination via SENSE
    img = SENSE_for_POCSENSE(ksp, RecMtx);
end
end