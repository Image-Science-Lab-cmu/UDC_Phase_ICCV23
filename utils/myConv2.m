function imgBlur = myConv2(img, kernels, shape)

if ~exist('shape', 'var')
    shape = 'same';
end

Isize = size(img);
ksize = size(kernels);
% sz = Isize - ksize + 1;

if length(size(img)) == 2
    numCh = 1;
else
    numCh = 3;
end

if numCh == 3
    % implement convolution
    % TODO: think about different boundary condition
    imgBlur = zeros(Isize(1)*2, Isize(2)*2, 3);
    for cc = 1:3
        kernel = kernels(:,:,cc);   % normalized done
        K = fft2(kernel, Isize(1)*2, Isize(2)*2);
        I = fft2(img(:,:,cc), Isize(1)*2, Isize(2)*2);
        imgBlur(:,:,cc) = ifft2(K.*I);
    end
    
    xx = floor(ksize(1)/2);
    yy = floor(ksize(2)/2);
    % todo
    switch shape
        case 'full'
            sz = Isize + ksize - 1;
            imgBlur = imgBlur(1: sz(1), 1: sz(2), :);
        case 'valid'
            sz = Isize - ksize + 1;
            imgBlur = imgBlur(2*xx: 2*xx+sz(1)-1, 2*yy:2*yy+sz(2)-1, :);
        otherwise
            imgBlur = imgBlur(xx+1:xx+Isize(1), yy+1:yy+Isize(2), :);
    end
    
else
    %% Gray-scale image
    
    % implement convolution
    max_sz = max(ksize, Isize);
    K = fft2(kernels, max_sz(1)*2, max_sz(2)*2);
    I = fft2(img, max_sz(1)*2, max_sz(2)*2);
    imgBlur = ifft2(K.*I);
    
    % IFFT2 when compute on Nvidia GPU
    % has decrepencies to CPU version.
    % There is a very small imagery component (~1e-17)
    % We manually ignore the imager component.
    if gpuDeviceCount > 0
        imgBlur = real(imgBlur);
    end
    
    xx = floor(ksize(1)/2);
    yy = floor(ksize(2)/2);
    % todo
    switch shape
        case 'full'
            sz = Isize + ksize - 1;
            imgBlur = imgBlur(1: sz(1), 1: sz(2));
        case 'valid'
            sz = Isize - ksize + 1;
            imgBlur = imgBlur(2*xx: 2*xx+sz(1)-1, 2*yy:2*yy+sz(2)-1);
        otherwise
            imgBlur = imgBlur(xx+1:xx+Isize(1), yy+1:yy+Isize(2));
    end
end