clc;
clear;
close all;
%more off;

% ------------------------------------------------------
% PART B EVALUATION
% ------------------------------------------------------

% --- INIT
if exist('OCTAVE_VERSION', 'builtin')>0
    % If in OCTAVE load the image package
    warning off;
    pkg load image;
    warning on;
end

% ------------------------------------------------------
% COMPARE RESULTS TO GROUND TRUTH
% ------------------------------------------------------

% --- Step B1
% load the ground truth
GT=readmatrix('Troizina 1827_ground_truth.txt');
% load our results (if available)
if exist('results.txt','file')
    R=readmatrix('results.txt');
else
    error('Ooooops! There is no results.txt file!');
end

% --- Step B2
% define the threshold for the IOU matrix
T=0.5; % 0.5 or 0.3 or 0.7

% calculate IOU for ;all the results
IOU = calcIOU(R,GT);

% apply the IOU threshold
IOUFinal = IOU >= T;


% calculate TP, FP, FN, Recall, Precision and F-Measure


bb_num = size(GT,1); % number of ground truth bounding boxes

% calculate TP 
TP = sum(IOUFinal(:));

% calculate FP 
FP = bb_num - TP;

% calculate FN
FN = bb_num - TP - FP;

% calculate recall
recall = TP / (TP + FN);

% calculate precision
precision = TP / (TP + FP);

% calculate f_measure
f_measure = 2 * (precision * recall) / (precision + recall);

% and show the results
fprintf('Recall %0.2f\n',recall);
fprintf('Precision %0.2f\n',precision);
fprintf('F-Measure %0.2f\n',f_measure);
