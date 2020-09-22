function [segimg,J,JF,Z,BW,BW2] = segmentBloodVessels(I)
%I = imread('1_org.tif');
% Resize image for easier computation
%B = imresize(I, [584 565]);
[r,c, ch] = size(I);
B = imresize(I, [r c]);
                    % Read image
gray = im2double(B);
%% Contrast Enhancment of gray image using CLAHE
J = adapthisteq(gray,'numTiles',[8 8],'nBins',128);
%% Background Exclusion                    % Apply Average Filter
h = fspecial('average', [11 11]);
JF = imfilter(J, h);
Z = imsubtract(JF, J);
%% Threshold using the IsoData Method
level = multithresh(Z);

 %% Convert to Binary
BW = im2bw(Z, level-.008);
%% Remove small pixels
 BW2 = bwareaopen(BW, 100);
 
%% Overlay
BW2 = imcomplement(BW2);
out = imoverlay(B, BW2, [0 0 0]);
           
vv = rgb2gray(out);
binaryImage = vv > 60;
                   
BW2 = bwmorph(binaryImage,'clean');
segimg = BW2;
end