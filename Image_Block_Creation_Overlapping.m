clc;
clear all;
Im = imread('a.jpg');
[M,N,~] = size(Im);
rr = 128; cc = 128; xx = 32; yy = 32;

numBlocksYY = numel(1:rr-xx:(M-(rr-1)));
numBlocksXX = numel(1:cc-yy:(N-(cc-1)));
C = cell(numBlocksYY*numBlocksXX,1);

counter = 1;
for ii=1:rr-xx:(M-(rr-1))
    for jj=1:cc-yy:(N-(cc-1))
        C{counter} =  Im(ii:(ii+rr-1), jj:(jj+cc-1), : );
        counter = counter + 1;
    end
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % indices of the first rr x cc block
% [YY,XX]=ndgrid( 1:(rr-xx):(M-(rr-1)) , 1:(cc-yy):(N-(cc-1)));
% IDX1 = bsxfun(@plus, (1:rr)', ((1:cc)-1)*M);
% % offset in indices for each block
% offset = (XX(:)-1) + (YY(:)-1)*M;
% % indices of each block, the block is indexed with the 3rd dimension
% IDX = bsxfun(@plus, IDX1, reshape(offset, [1 1 numel(offset)]));
% % convert to cell of blocks
% C = mat2cell(Im(IDX), rr, cc, ones(1, numel(XX)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure;
for ii=1:numBlocksYY*numBlocksXX
    subplot(numBlocksYY,numBlocksYY,ii), imshow( C{ii} ); axis image; colormap gray;
end
%%
for ii=1:numBlocksYY*numBlocksXX
    nnzCount(ii) = nnz(C{ii});
end

[ma,loc]=max(nnzCount)
