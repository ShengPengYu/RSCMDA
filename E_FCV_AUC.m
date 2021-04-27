%Prediction work
clc;              
clear;            % clear all workspace variables
close all;        % close all windows
  
% matlab does not have the concept of layering, so the data from other 
% subfolders is added to the main program directory before adding code
% to the program:addpath(genpath(pwd));
currentFolder = pwd;              
addpath(genpath(currentFolder));   
load SD;
AD = SD;
load SM;
AM = SM;

load knownre ;
Y = knownre ;         % Y is the ground truth matrix (383*495)
load HMDD;
load x ;
X = x ;


global_position =importdata('./Experiments/fcv/fcv_position.txt');
F_AUC = Fpositiontooverallauc(Y,HMDD,global_position');

