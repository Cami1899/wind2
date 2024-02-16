clear all
close all
clc
format long

pathstr = fileparts(mfilename('fullpath'));
addpath(sprintf('%s/Functions',pathstr));
addpath(sprintf('%s/OdeFunction',pathstr));
addpath(sprintf('%s/Script',pathstr));
addpath(sprintf('%s/Data',pathstr));
addpath(sprintf('%s/WindModel',pathstr));
savepath;

cd(pathstr);


load('wind_value_restricted.mat');
load('time_restricted.mat');