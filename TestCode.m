clc;
clear all;
tic;
I = imread('Image Dataset\test\images\01_test.tif');
figure; imshow(I);
grayI = rgb2gray(I);
grayI = imsharpen(grayI);
BW = imbinarize(grayI);
grayI = grayI.*uint8(BW);
imwrite(grayI,'Output Images\01_BW.jpg');

[segmentedBloodVessel,J,JF,Z,BW,BW2] = segmentBloodVessels(grayI);
% imwrite(J,'Output Images\01_CLAHE.jpg');
% imwrite(JF,'Output Images\01_Filter.jpg');
% imwrite(Z,'Output Images\01_BackgroundExclusion.jpg');
% imwrite(BW,'Output Images\01_Binary.jpg');
% imwrite(BW2,'Output Images\01_NoiseRemoved.jpg');

%imtool(segmentedBloodVessel);
%imwrite(I,'Output Images\01_Original.jpg');
%imwrite(segmentedBloodVessel,'Output Images\01_Blood_Vessel_Mask.jpg');


%%
Im = segmentedBloodVessel;
[M,N,~] = size(Im);
rr = 132; cc = 132; xx = 32; yy = 32;

numBlocksYY = numel(1:rr-xx:(M-(rr-1)));
numBlocksXX = numel(1:cc-yy:(N-(cc-1)));
C = cell(numBlocksYY*numBlocksXX,1);
CLoc = cell(numBlocksYY*numBlocksXX,1);

counter = 1;
for ii=1:rr-xx:(M-(rr-1))
    for jj=1:cc-yy:(N-(cc-1))
        C{counter} =  Im(ii:(ii+rr-1), jj:(jj+cc-1), : );
        CLoc{counter} = [ii jj];
        counter = counter + 1;
    end
    fprintf('\n');
end

%%
figure;
for ii=1:numBlocksYY*numBlocksXX
    subplot(numBlocksYY,numBlocksXX,ii), imshow( C{ii} ); axis image; colormap gray;
end

%%
nnzCount = zeros(numBlocksYY*numBlocksXX,1);
for ii=1:numBlocksYY*numBlocksXX
    nnzCount(ii,1) = nnz(C{ii});
end

[ma,loc]=max(nnzCount);

%%
pos = CLoc{loc};
CropI = I(pos(1):(pos(1)+rr-1), pos(2):(pos(2)+cc-1), : );
figure;imshow(CropI)
%imwrite(CropI,'Output Images\01_CropI.jpg');

%%
f = fspecial('average', [5 5]);
filtCropI = imfilter(CropI(:,:,1),f);
Z = imsubtract(CropI(:,:,1),imcomplement(filtCropI));
level = graythresh(Z);
mask = imbinarize(Z,level);
SE = strel('disk',25);
mask = imopen(mask,SE);
figure;imshow(mask) 

%%
NewMask = zeros(M,N);
NewMask(pos(1):(pos(1)+rr-1), pos(2):(pos(2)+cc-1), : ) = mask;
imshow(NewMask)
%imwrite(NewMask,'Output Images\01_ODmask.jpg');

%%
B = bwboundaries(NewMask);
stat = regionprops(NewMask,'Centroid');
figure;imshow(I); hold on
for k = 1 : length(B)
    b = B{k};
    c = stat(k).Centroid;
    %c(1) = c(1)+40; 
    viscircles(c,40,'Color','g');
    text(c(1),c(2),num2str(k),'backgroundcolor','g');
end
toc;