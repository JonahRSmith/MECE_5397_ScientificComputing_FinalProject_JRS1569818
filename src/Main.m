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
DTIMEI=0.0125; %Explicit routine requires smaller DTIMEI than ADI for convergence to steady state
maxrelerror=1e-14; %Max relative error to classify steady state
savefilename='Test_Save.mat';
%1/DTIMEI = Frequency a save is made, and that convergence is checked
%Call solution routine
if strcmp(ui,'ADI')
    [x,y,u,runtime]=ADI(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);
else
    [x,y,u,runtime]=Explicit(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);
end
if max(max(isnan(u)))
    warning('Solution for u has NaN values in it. Reduce DTIMEI to converge to a real solution.');
end
surf(x,y,u');
% %% Calling ADI using a structure
% struct_file = load('SaveFile.mat');
% [x,y,u]=ADI(struct_file);

% %% Running a loop around Explicit.m and ADI.m
% clc; clear all;
% ax=0; bx=2*pi;
% ay=0; by=2*pi;
% DTIMEI=0.04;
% maxrelerror=1e-14; %Max relative error to classify steady state
% for refine=3:7
%     fprintf('In refine = %i\n',refine);
%     DTIMEI=DTIMEI/4
%     nodefacx=refine;
%     nodefacy=refine;
%     savefilename=['Save_Explicit_',num2str(refine),'.mat'];
%     [x,y,u,runtime]=Explicit(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);
%     if max(max(isnan(u)))
%         warning('Solution for u has NaN values in it. Reduce DTIMEI to converge to a real solution.');
%     end
% end
% 
% ax=0; bx=2*pi;
% ay=0; by=2*pi;
% DTIMEI=0.4;
% maxrelerror=1e-14; %Max relative error to classify steady state
% for refine=3:7
%     fprintf('In refine = %i\n',refine);
%     DTIMEI=DTIMEI/4
%     nodefacx=refine;
%     nodefacy=refine;
%     savefilename=['Save_ADI_',num2str(refine),'.mat'];
%     [x,y,u,runtime]=Explicit(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);
%     if max(max(isnan(u)))
%         warning('Solution for u has NaN values in it. Reduce DTIMEI to converge to a real solution.');
%     end
% end


%% Plotting
surf(x,y,u');
colorbar;
view(0,90);
% Report DTIMEI and relerror for some run
% load('YourSaveFile.mat');
% DTIMEI
% xnodes=size(x); ynodes=size(y);
% relerror=0;
% for j=1:ynodes
%     for i=1:xnodes
%         relerror=max(relerror,abs((u(j,i)-uprev(j,i))/u(j,i)));
%     end
% end
% relerror

%% Comparing two saves
S1 = load('Save_Explicit_3.mat');
u1=S1.u;
clc;

[xn,yn]=size(u1);
max1=max(max(u1(2:xn,2:yn-1)));
min1=min(min(u1(2:xn,2:yn-1)));
mean1=mean(mean(u1(2:xn,2:yn-1)));
fprintf('%g\t%g\t%g\n',max1,min1,mean1);

%% Explicit_Verification.m
%Verification using a simple 4-Dirichlet boundary u
% u(ax,y) = 0;
% u(bx,y) = 1;
% u(x,ay) = (x-ax)/(bx-ax);
% u(x,by) = (x-ax)/(bx-ax);
%We should see a flat plane tilted at an angle with a slope 1/(bx-ax)
ax=0; bx=2*pi;
ay=0; by=2*pi;
nodefacx=6;
nodefacy=6;
DTIMEI=0.00125; %Explicit routine requires smaller DTIMEI than ADI for convergence to steady state
maxrelerror=1e-14; %Max relative error to classify steady state
savefilename='Test_Save.mat';
[x,y,u,runtime]=Explicit_Verification(ax,bx,ay,by,nodefacx,nodefacy,DTIMEI,maxrelerror,savefilename);

surf(x,y,u');




