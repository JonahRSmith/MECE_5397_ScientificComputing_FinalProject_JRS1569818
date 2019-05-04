function [ x,y,u ] = ADI( varargin )
%Solve a 2D diffusion equation using the ADI method
%Solution is ran to steady state
%Saves of the workspace are made at regular intervals
% Proper number of inputs provided to function
%   varargin(1) = ax
%	varargin(2) = bx
%	varargin(3) = ay
%	varargin(4) = by
%	varargin(5) = nodefacx ~ Defines level of nodalization of problem in x
%	              Increasing nodefac means exponentially more number of
%	              nodes in the x-axis
%	varargin(6) = nodefacy ~ Defines level of nodalization of problem in y
%	              Increasing nodefac means exponentially more number of
%	              nodes y-axis
%	varargin(7) = DTIMEI ~ Defines the timestep taken each iteration by
%	              the solution routine. Decreasig DTIMEI decreases
%	              the amount of time elapsed each iteration. Setting
%	              DTIME <= 0 will cause the program to automatically
%	              calculate a value of DTIMEI for you.
%   varargin(8) = maxrelerror ~ Defines the maximum relative error in
%                 any entry in u that is allowed between two iterations
%                 to classify convergence to steady state as having been
%                 reached.
%   varargin(9) = savefilename ~ Name that is saved to.
%By passing through a saved structure from a prevous run as the only
%argument to this function, the solution will start from the end of the
%save and results will be appended to the previously ran solution.


desired_inputs = 9;
tic;
fprintf('Setting up problem to perform ADI routine.\n');

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
    xnodes = (ceil(4^(nodefacx)));
    ynodes = (ceil(4^(nodefacy)));
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
    fb = (bx-x).^2 .* cos(pi*x/bx); %Vector containing fb(x) at each x
    gb = x.*(bx-x).^2; %Vector containing gb(x) at each x
    for xi=1:xnodes %Top and bottom Dirichlet conditions
        u(xi,1) = gb(xi);
        u(xi,ynodes) = fb(xi);
        %Initial guess: linear relationship between top and bottom
        for yi=2:ynodes-1
            u(xi,yi) = ((ay+DY*(yi-1))/(ay+by))*(u(xi,ynodes)-u(xi,1))+u(xi,1);
        end
        %Though it might be better just to assume a 0
    end
    %Dirichlet condition on x=ax boundary
    for yi=1:ynodes
        %gb(ax) + (y-ay)/(by-ay)*(fb(ax)-gb(ax))
        u(1,yi)=gb(1) + (y(yi)-ay)/(by-ay)*(fb(1)-gb(1));
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
uprev=u;
for j=2:ynodes-1
    for i=2:xnodes-1
        u(i,j) = uprev(i,j) + DTIMEIperDX2*(uprev(i-1,j)-2*uprev(i,j)+uprev(i+1,j)) + DTIMEIperDY2*(uprev(i,j-1)-2*uprev(i,j)+uprev(i,j+1));
    end
    %Now solving at i=xnodes, the Neumann boundary
    u(xnodes,j) = uprev(xnodes,j) + DTIMEIperDX2*(2*uprev(xnodes-1,j)-2*uprev(xnodes,j)) + DTIMEIperDY2*(uprev(xnodes,j-1)-2*uprev(xnodes,j)+uprev(xnodes,j+1));
end
    
%=-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-==-=-=-=-=-=
    TIMEN = TIMEN+DTIMEI;
    %fprintf('TIMEN=%g\n',TIMEN);
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
end
fprintf('Convegence met at TIMEN=%g\n',TIMEN);














end