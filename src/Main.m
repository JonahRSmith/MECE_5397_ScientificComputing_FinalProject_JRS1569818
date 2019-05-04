%% Project B - Diffusion Equation
%Jonah R. Smith, 1569818

%This file serves as a wrapper to setup the problem
%The problem is solved in one of the following functions:
%Explicit.m, or ADI.m
%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
clc; clear all; close all;
%Load in parameters

ax=0; bx=2*pi;
ay=0; by=2*pi;
nodefacx=3;
nodefacy=3;
DTIMEI=0.01;
maxrelerror=1.2e-15; %Max relative error to classify steady state
savefilename='Save001.mat';
%1/DTIMEI = Frequency a save is made, and that convergence is checked
%Call solution routine
[x,y,u]=ADI(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);
% %% Calling ADI using a structure
% struct_file = load('SaveFile.mat');
% [x,y,u]=ADI(struct_file);