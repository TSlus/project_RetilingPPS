clear;clc;
addpath(genpath('Models'));

name1 = 'polyhedron.mat'; 
name2 = 'fandisk.mat'; 
name3 = 'block.mat'; 
name4 = 'star.mat'; 
%% input
% 1.mesh
modelname = name1;
load(modelname);%obj_write('block.obj', vertices, faces);
modelname = 0;% when you want to use the saved data, need it.

% 2.number of candidate points
nCand = ceil(size(vertices, 1));

% 3.the level of detail, high curvature old points -- CSP
% k_level = 0;
k_level = ceil(nCand*3/5); % quarter of nCand

% 4.methods for remove old vertices
method = 3; % 1,2,3，both = 1,2,3，4

%% do Retiling and PPS
detail_plot = 0;
Retiling_and_PPS;
