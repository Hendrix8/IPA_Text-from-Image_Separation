clc;
clear;
close all;

% ------------------------------------------------------
%                   PART A MAIN 
% ------------------------------------------------------

% --- INIT
if exist('OCTAVE_VERSION', 'builtin')>0
    % If in OCTAVE load the image package
    warning off;
    pkg load image;
    warning on;
end

% ------------------------------------------------------
%   LOAD AND SHOW THE IMAGE + CONVERT TO BLACK-AND-WHITE
% ------------------------------------------------------

% --- Step A1
% read the original RGB image 
Filename='Troizina 1827.jpg';
I=imread(Filename);

% show it (Eikona 1)
figure;
image(I);
axis image off;
colormap(gray(256));

% --- Step A2
% convert the image to grayscale
A=any_image_to_grayscale_func('Troizina 1827.jpg');

% apply gamma correction (a value of 1.0 doesn't change the image)
GammaValue=1.0; 
A=imadjust(A,[],[],GammaValue); 

% show the grayscale image (Eikona 2)
figure;
image(A);
colormap(gray(256));
axis image off;
title('Grayscale image');

% --- Step A3
% convert the grayscale image to black-and-white using Otsu's method
Threshold= graythresh(A); 
BW = imbinarize(A,Threshold);
imwrite(BW, 'images/bw.png');
% show the black-and-white image (Eikona 3)
figure;
image(BW);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);
title('Binary image');

% ------------------------------------------------------
%CLEAN THE IMAGE
% ------------------------------------------------------

% --- Step A4
% making morphological operations to clean the image 

% ------------------------------------------------------
%                   DILATING THE IMAGE
% ------------------------------------------------------


se_dil = strel('disk', 4); % creating a structuring element 
clean_im = ~imdilate(~BW, se_dil); % dilating the image using the dilation structuring element

% showing the dilated image
figure;
image(clean_im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/dil.png');

% ------------------------------------------------------------------------------------------
%           FINDING THE BIGGEST CONNECTED COMPONENT AND REMOVING IT
% ------------------------------------------------------------------------------------------

clean_im = ~remove_bigComp(~clean_im);

% showing the image with the removed stamp
figure;
image(clean_im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/stmpRem.png');


% ------------------------------------------------------------------------------------------
%                     ADDING PADDING AND THEN REMOVING THE BORDER WITH 
%                    THE SAME LOGIC AS IT WILL THEN BE THE BIGGEST COMPONENT
% ------------------------------------------------------------------------------------------

clean_im = padarray(clean_im, [20 20], 0); % adding padding to the image
imwrite(clean_im, 'images/pad.png');

se_dil = strel('disk', 2); % creating a structuring element 
clean_im = ~imdilate(~clean_im, se_dil); % dilating the image a little bit for the border to connect

 % showing the image with the padding
figure;
image(clean_im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/er_enhance_pad.png');

% removing the border with the function remove_bigComp
clean_im = ~remove_bigComp(~clean_im);

[rows, cols, ~] = size(clean_im);
clean_im = imcrop(clean_im, [20 20 cols-40 rows-40]);

% showing the image with the removed padding
figure;
image(clean_im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/bord_rem.png')

% ------------------------------------------------------------------------------------------
%                          ERODING THE IMAGE TO REMOVE THE NOISE 
% ------------------------------------------------------------------------------------------

se_er = strel('disk', 8); % creating a structuring element 
clean_im = ~imerode(~clean_im, se_er); % eroding the image using the eroding structuring element

% showing the image after the eroding
figure;
image(clean_im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/noise_rem.png')

% ------------------------------------------------------------------------------------------
%                          DILATING THE IMAGE TO ENHANCE THE WORDS
% ------------------------------------------------------------------------------------------

se_dil = strel('disk', 3); % creating a structuring element 
clean_im = ~imdilate(~clean_im, se_dil); % dilating the image using the dilating structuring element

 % showing the image after the dilating
figure;
image(clean_im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/enhance_words.png')

% ------------------------------------------------------------------------------------------
%                       REMOVING THE PART THAT IS IN BETWEEN ΤΗΕ TITLE AND
%                       SUBTITLE
% ------------------------------------------------------------------------------------------

[rows, cols] = size(clean_im); % getting the metrics of the image

mask = ones(rows, cols); % creating a mask that has ones and is the same size as the image

% making the mask so that it removes the part in between
mask(round(rows/6):round(rows/4.7),:) = 0; 

clean_im = immultiply(~clean_im,mask); % multiplying the original image with the mask so that the every part except the upper half is erased

figure;
image(~clean_im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);
print(gcf, '-r300', '-dpng', 'images/removed_part.png')


% ------------------------------------------------------
% WORD SEGMENTATION
% ------------------------------------------------------
% --- Step A5

% make morphological operations for word segmentation ...

title_mask = zeros(rows, cols); % creating a mask that has zeros and is the same size as the image
title_mask(1:round(0.5*rows/3),:) = 1; % cutting the mask to show only the title

title = immultiply(clean_im,title_mask); % multiplying the original image with the mask so that the every part except the upper half is erased

figure;
image(~title);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng','images/onlytitle.png')

subtitle_mask = zeros(rows, cols); % mask of zeros
subtitle_mask(round(0.5*rows/3):round(rows/3),:) = 1; % replacing with ones the subtitle part

subtitle = immultiply(clean_im, subtitle_mask); % applying the mask

figure;
image(~subtitle);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng','images/onlysub.png');

subtitle2_mask = zeros(rows, cols); % mask of zeros
subtitle2_mask(round(rows/3):380,:) = 1; % replacing with ones the subtitle part

subtitle2 = immultiply(clean_im, subtitle2_mask); % applying the mask

figure;
image(~subtitle2);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng','images/onlysub2.png');

% dilating the cropped images
se_dil = strel('line',25, 0); % creating a structuring element of line 
title = imdilate(title, se_dil); % dilating the title

se_dil = strel('line', 15, 0); % creating a structuring element of line 
subtitle = imdilate(subtitle, se_dil); % dilating the subtitle

se_dil = strel('line', 28, 0); % creating a structuring element of line
subtitle2 = imdilate(subtitle2, se_dil); % dilating the subtitle

figure;
image(~title);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/only_title_dilated.png')


figure;
image(~subtitle);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);
print(gcf, '-r300', '-dpng', 'images/subonly_title_dilated.png')

figure;
image(~subtitle2);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);
print(gcf, '-r300', '-dpng', 'images/subonly2_title_dilated.png')
 
% adding cropped images to the original image
im = clean_im | title; % adding the dilated title
im = im | subtitle; % adding the dilated subtitle
im = im | subtitle2; % adding the last subtitle

figure;
image(~im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);

print(gcf, '-r300', '-dpng', 'images/title_dilated.png');


im(445:455,:) = 0; % cutting the connections between lines in a certain spot

figure;
image(~im);
colormap(gray(2));
axis image;
set(gca,'xtick',[],'ytick',[]);
print(gcf, '-r300', '-dpng', 'images/image_seg.png');
 
% find the connected components (using bwlabel) 

labels = bwlabel(im, 4); % finding all the connected components with connectivity=4
label4_num = max(labels(:)) % number of labels found

labels8 = bwlabel(im, 8);
label8_num = max(labels(:))

figure;
image(label2rgb(labels, 'lines'))
axis image;
set(gca,'xtick',[],'ytick',[]);
print(gcf, '-r300', '-dpng', 'images/labels.png');

% figure;
% image(label2rgb(labels8, 'lines'))
% axis image;
% set(gca,'xtick',[],'ytick',[]);
% print(gcf, '-r300', '-dpng', 'images/labels8.png');

% ------------------------------------------------------
% FINAL IMAGE WITH BOUNDING BOXES
% ------------------------------------------------------
% --- Step A6


% --- Store all the bounding boxes in an array R
R=[];
for i=1:label4_num
    [r,c]=find(labels==i);
    XY=[c r];
    x1=min(XY(:,1));
    y1=min(XY(:,2));
    x2=max(XY(:,1));
    y2=max(XY(:,2));
    R=[R;x1 y1 x2 y2]; % append to R the bounding box as [x1 y1 x2 y2]
end

% -------------------------------------------------------------------------------------------
%                   ALTERING THE BOUNDING BOXES TO INCLUDE PUNCTUATION
%                   (Comment to see original bounding boxes)
% -------------------------------------------------------------------------------------------

R=[];
for i=1:label4_num
    [r,c]=find(labels==i);
    XY=[c r];
    x1=min(XY(:,1)) - 8;
    y1=min(XY(:,2)) - 8;
    x2=max(XY(:,1)) + 3;
    y2=max(XY(:,2));
    R=[R;x1 y1 x2 y2]; 
end



% show word segmentation (Eikona 5)

% --- Draw the original image with the bounding boxes around the connected
% --- components
ColorMap=lines;
figure('color','w');
A=(I);
image(A); % NOTICE! the image is shown inverted
colormap(gray(256));
axis image;
set(gca,'xtick',[],'ytick',[]);
title('Bounding boxes');


for i=1:size(R,1)
    x1=R(i,1);
    y1=R(i,2);
    x2=R(i,3);
    y2=R(i,4);
    x = [x1-0.5 x2+0.5 x2+0.5 x1-0.5 x1-0.5];
    y = [y1-0.5 y1-0.5 y2+0.5 y2+0.5 y1-0.5];
    line(x, y, 'color',ColorMap(mod(i-1,7)+1,:),'linewidth',0.5);
end

print(gcf, '-r300', '-dpng', 'images/final_improved.png');






% --- Step A7
% let suppose matrix R contains all the N bounding boxes in the form
% x11 y11 x12 y12
% x21 y21 x22 y22
% ...
% xN1 yN1 xN2 yN2
% then save the bounding boxes in a text file results.txt ...

%s=regionprops(labels,'BoundingBox'); % properties of the bounding boxes
%AllBoxes=cat(1,s.BoundingBox); % saving the metrics of every bounding box in a matrix
if exist('R','var')
     dlmwrite('results.txt',R,'\t');
 else
     error('Ooooops! There is no R variable!');
end


%----------------------------------------------------------
%                       FUNCTIONS
%----------------------------------------------------------

% this function takes the components of an image and returns the index of
% the index of the biggest component
function m = max_comp(cc_areas)
    max_area = 0;
    max_component = 0;
    for i=1:length(cc_areas)
      if cc_areas(i).Area > max_area
        max_area = cc_areas(i).Area;
        max_component = i;
      end
    end
    m = max_component;
end

% Given a binary image, this function removes the
% biggest connected component 
function I = remove_bigComp(bw)

    conn_comps = bwconncomp(bw, 8); % finding all the connected components with connectivity = 4
    cc_areas = regionprops(conn_comps, 'Area'); % finding the area of all the connected components of the image
    
    % finding the index of the biggest component using the max_comp function
    max_idx = max_comp(cc_areas);
    
    % removing the biggest component from the image
    bw(conn_comps.PixelIdxList{max_idx}) = 0;
    I = bw;

end



