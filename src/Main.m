%% Project B - Diffusion Equation
%Jonah R. Smith, 1569818

%This file serves as a wrapper to setup the problem
%The problem is solved in one of the following functions:
%Explicit.m, or ADI.m
%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
clc; clear all; close all;
%Load in parameters
ui=questdlg('Which solution method would you like to use?','Choose Solution Routine','ADI','Explicit','Explicit');

ax=0; bx=2*pi;
ay=0; by=2*pi;
nodefacx=6;
nodefacy=6;
DTIMEI=0.01; %Explicit routine requires smaller DTIMEI than ADI for convergence to steady state
maxrelerror=1e-14; %Max relative error to classify steady state
savefilename='Save_Explicit_001.mat';
%1/DTIMEI = Frequency a save is made, and that convergence is checked
%Call solution routine
if strcmp(ui,'ADI')
    [x,y,u]=ADI(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);
else
    [x,y,u]=Explicit(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);
end
if max(max(isnan(u)))
    warning('Solution for u has NaN values in it. Reduce DTIMEI to converge to a real solution.');
end
% %% Calling ADI using a structure
% struct_file = load('SaveFile.mat');
% [x,y,u]=ADI(struct_file);