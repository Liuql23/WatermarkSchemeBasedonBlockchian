function [cost]=DistortionCost(cover)

hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, ...
        -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768]; % 1D high pass decomposition filter
lpdf = (-1).^(0:numel(hpdf)-1).*fliplr(hpdf);% 1D low pass decomposition filter
% construction of 2D wavelet filters
F{1} = lpdf'*hpdf;
F{2} = hpdf'*lpdf;
F{3} = hpdf'*hpdf;


padSize = max([size(F{1})'; size(F{2})'; size(F{3})']);
sgm=1;
[lm,ln]=size(cover);
coverPadded = padarray(cover, [padSize padSize], 'symmetric');
xi = cell(3, 1);
for fIndex = 1:3
    % compute residual
    R = conv2(coverPadded, F{fIndex}, 'same');
    % compute suitability
    xi{fIndex} = conv2(1./(abs(R)+sgm), rot90(abs(F{fIndex}), 2), 'same');
    % correct the suitability shift if filter size is even
    if mod(size(F{fIndex}, 1), 2) == 0, xi{fIndex} = circshift(xi{fIndex}, [1, 0]); end;
    if mod(size(F{fIndex}, 2), 2) == 0, xi{fIndex} = circshift(xi{fIndex}, [0, 1]); end;
    % remove padding
    xi{fIndex} = xi{fIndex}(((size(xi{fIndex}, 1)-lm)/2)+1:end-((size(xi{fIndex}, 1)-lm)/2), ((size(xi{fIndex}, 2)-ln)/2)+1:end-((size(xi{fIndex}, 2)-ln)/2));
end
cost= xi{1} + xi{2} + xi{3};
end