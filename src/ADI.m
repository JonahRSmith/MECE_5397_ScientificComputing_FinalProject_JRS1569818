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
%By passing through a saved structure from a prevous run as the only
%argument to this function, the solution will start from the end of the
%save and results will be appended to the previously ran solution.


desired_inputs = 7;

fprintf('Setting up problem to perform ADI routine.\n');

if nargin < 1, error('Must provide arguments when calling ADI.m'); end
if nargin == 1,
    if ~isstruct(varargin(1))
        error('Save files must be provided as a structure containing all variables when calling CrankNicolson.n');
    else
        %Load in savefile from previous run
        varlist = fieldnames(varargin(1));
        for vn=1:length(varlist)
            cmdstr=[varlist(vn),'=varlist(vn);'];
            eval(cmdstr); %Unpacks each variable into the workspace
        end
        clear varlist vn cmdstr;
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
    %Steady state solution time would be quicker if I found an initial
    %condition closer to the steady state solution, but the purpose of this
    %code is to be able to find that steady state solution regardless of
    %the initial condition, so this approximation will do.
    %Tri-Diagonal Solver (ADI)
    %   Recall: We add one column of ghost nodes over x at x=bx+DX
    for zz=1:750
        %Iterate over x with half a timestep then over y
        uprev = u; %Solution of previous timestep
        N=xnodes-2+1; %1 for the ghost node
        for j=2:ynodes-1
            alpha=zeros(N);
            g=zeros(N);
            alpha(1)=(2/DTIMEI + 2/DX^2);
            g(1)=(2/DTIMEI - 2/DY^2)*uprev(2,j) + (1/DY^2)*(uprev(2,j-1)+uprev(2,j+1));
            norma=alpha(1);
            normb = -1/DX^2; %Saving this so I don't redo it over and over
            normb2 = normb^2;
            for i=2:N-1 % -1 so I can fix the ghost node
                alpha(i) = norma - normb2/alpha(i-1);
                g(i)=(2/DTIMEI - 2/DY^2)*uprev(i+1,j) + (1/DY^2)*(uprev(i+1,j-1)+uprev(i+1,j+1)) - normb*g(i-1)/alpha(i-1);
            end
            %Ghost node:
            alpha(N) =  norma - 2*normb2/alpha(N-1);
            g(N)=(2/DTIMEI - 2/DY^2)*uprev(N+1,j) + (1/DY^2)*(uprev(N+1,j-1)+uprev(N+1,j+1)) - 2*normb*g(N-1)/alpha(N-1);
            u(N+1,j) = g(N)/alpha(N);
            for k=N-1:-1:1
                u(k+1,j) = (g(k)-normb*u(k+2,j))/alpha(k);
            end
        end
        %Now to iterate over y
        uprev=u;
        N=ynodes-2;
        for j=2:xnodes-1
            alpha=zeros(N);
            g=zeros(N);
            alpha(1)=(2/DTIMEI + 2/DY^2);
            g(1)=(2/DTIMEI - 2/DY^2)*uprev(j,2) + (1/DY^2)*(uprev(j-1,2)+uprev(j+1,2));
            norma=alpha(1);
            normb = -1/DY^2; %Saving this so I don't redo it over and over
            normb2 = normb^2;
            for i=2:N
                alpha(i) = norma - normb2/alpha(i-1);
                g(i)=(2/DTIMEI - 2/DY^2)*uprev(j,i+1) + (1/DY^2)*(uprev(j-1,i+1)+uprev(j+1,i+1)) - normb*g(i-1)/alpha(i-1);
            end
            u(j,N+1) = g(N)/alpha(N);
            for k=N-1:-1:1
                u(j,k+1) = (g(k)-normb*u(j,k+2))/alpha(k);
            end
        end
        fprintf('zz=%i\n',zz);
    end
end
















end