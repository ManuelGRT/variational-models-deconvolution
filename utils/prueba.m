clear all, close all, clc

basePath = fileparts(pwd);
img_folder = fullfile(basePath, 'images');
im = imread(fullfile(img_folder, 'panda.png'));
im = im2double(im); % im2double transforms the data into double type and in the range [0,1]

imshow(im)
