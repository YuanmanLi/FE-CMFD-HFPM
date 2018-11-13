
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% This code was developed by Li Yuanman 
% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite the paper if you use this code:
% Y. M. Li and J. T. Zhou, ¡°Image Copy-Move Forgery Detection Using Hierarchical Feature Point Matching¡±, ASC, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc 
clear all
dbstop if error
addpath 'lib\vlfeat-0.9.18\mexw64' 'functions'

img_name = '00001_2_scale.jpg';

im_file = ['test_pic\', img_name];
fprintf('Processing: %s\n',im_file);
tstart = tic;
%%input: 
%%first para: input file
%%second para: whether show the matches, default 0
%%third para: whether conduct large resizing (set 1 when s>2). default 0. 
[countTrasfGeom,map] = process_image(im_file,0,1);
tproc = toc(tstart);
tps = datestr(datenum(0,0,0,0,0,tproc),'HH:MM:SS');

%imwrite(map, ['test_pic\result_', img_name])
imshow(map)
fprintf('\nComputational time: %s\n', tps);