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
nodefacx=4;
nodefacy=4;
DTIMEI=0.00001;
%Call solution routine
[x,y,u]=ADI(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI);