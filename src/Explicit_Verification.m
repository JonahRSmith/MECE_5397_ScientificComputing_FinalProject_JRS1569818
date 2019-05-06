function [ x,y,u,runtime ] = Explicit_Verification( varargin )
% This is a solution to a four-sided Dirichlet boundary Diffusion equation
% used for verification purposes. This IS NOT a solution to the given
% poblem, but instead, a solution to a known problem. The Explicit Method
% was used here rather than ADI as it is easy to set up an additional
% Dirichlet boundary in this script.
% Definition of Boundaries:
% u(ax,y) = 0;
% u(bx,y) = 1;
% u(x,ay) = (x-ax)/(bx-ax);
% u(x,by) = (x-ax)/(bx-ax);

desired_inputs = 9;
tic;
fprintf('Setting up problem to perform Explicit_Verification routine.\n');

if nargin < 1, error('Must provide arguments when calling ADI.m'); end
if nargin == 1,
    if ~isstruct(varargin{1})
        error('Save files must be provided as a structure containing all variables when calling ADI.m');
    else
        %Load in savefile from previous run
        SaveStruct=varargin{1};
        varlist = fieldnames(SaveStruct);
        for vn=1:length(varlist)
            cmdstr=[varlist{vn},'=SaveStruct.',varlist{vn},';'];
            eval(cmdstr); %Unpacks each variable into the workspace
        end
        clear varlist vn cmdstr SaveStruct;
    end
elseif nargin ~= desired_inputs
    error(sprintf('Must provide %i inputs, or one structure containing previous run''s information.',desired_inputs));
else
    %Loading in base variables
    ax = varargin{1};
    bx = varargin{2};
    ay = varargin{3};
    by = varargin{4};
    nodefacx = varargin{5};
    nodefacy = varargin{6};
    DTIMEI = varargin{7};
    maxrelerror = varargin{8};
    savefilename = varargin{9};
    
    %Calculate number of internal nodes in X and Y using the following formula
    xnodes = (ceil(2^(nodefacx)));
    ynodes = (ceil(2^(nodefacy)));
    DX = (bx-ax)/(xnodes-1);
    DY = (by-ay)/(ynodes-1);
    x = zeros(xnodes,1);
    y = zeros(ynodes,1);
    for xi=0:xnodes-1
        x(xi+1) = ax+DX*xi;
    end
    for yi=0:ynodes-1
        y(yi+1) = ay+DY*yi;
    end
    u = zeros(xnodes,ynodes);
    
    %Setting up problem
    TIMEN = 0; %Current time of the problem
    %Setting up guess for initial value of u
    fb = (x-ax)/(bx-ax); %Vector containing fb(x) at each x
    gb = (x-ax)/(bx-ax); %Vector containing gb(x) at each x
    %Setting up y=ax and y=bx Dirichlet conditions
    for xi=1:xnodes %Top and bottom Dirichlet conditions
        u(xi,1) = gb(xi);
        u(xi,ynodes) = fb(xi);
        %Though it might be better just to assume a 0
    end
    %Setting up Dirichlet condition on x=ax boundary and x=bx
    for yi=1:ynodes
        u(1,yi)=0;
        u(xnodes,yi)=1;
    end
    clear yi
end
%Steady state solution time would be quicker if I found an initial
%condition closer to the steady state solution, but the purpose of this
%code is to be able to find that steady state solution regardless of
%the initial condition, so this approximation will do.
%Tri-Diagonal Solver (ADI)
%   Recall: We add one column of ghost nodes over x at x=bx+DX
savefreq=ceil(1/DTIMEI)*50;
convergence=0;
fprintf('Beginning Solution Routine:\n');
DTIMEIperDX2 = DTIMEI/DX^2;
DTIMEIperDY2 = DTIMEI/DY^2;
while ~convergence
for zz=1:savefreq %Solving this many iterations before savings and checking for convergence
%=-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-=
%Explicit solver is extremely simple, as each node is solved independently
%of the others within a single timestep. This does however require a more
%refined DTIMEI in order to not have issues with divergence or oscillations
%in failure to converge.
uprev=u; %uprev holds the solution for u from the previous timestep
for j=2:ynodes-1 %Because I set up u to contain the boundaries inside of it,
    for i=2:xnodes-1 %No extra code is needed to apply Dirichlet conditions
        u(i,j) = uprev(i,j) + DTIMEIperDX2*(uprev(i-1,j)-2*uprev(i,j)+uprev(i+1,j)) + DTIMEIperDY2*(uprev(i,j-1)-2*uprev(i,j)+uprev(i,j+1));
    end
    %Now solving at i=xnodes, the Neumann boundary
    % u(xnodes,j) = uprev(xnodes,j) + 2*DTIMEIperDX2*(uprev(xnodes-1,j)-uprev(xnodes,j)) + DTIMEIperDY2*(uprev(xnodes,j-1)-2*uprev(xnodes,j)+uprev(xnodes,j+1));
end
TIMEN = TIMEN+DTIMEI;
%=-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-=
end
%Check for convergence
relerror=0;
for j=1:ynodes
    for i=1:xnodes
        relerror=max(relerror,abs((u(j,i)-uprev(j,i))/u(j,i)));
    end
end
fprintf('TIMEN=%g; Max relative error is %g; ',TIMEN,relerror); toc;
convergence=(relerror<=maxrelerror);
save(savefilename);
runtime=toc;
end
fprintf('Convegence met at TIMEN=%g\n',TIMEN);














end