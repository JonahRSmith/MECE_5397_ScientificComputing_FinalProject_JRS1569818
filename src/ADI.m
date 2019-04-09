function [ output_args ] = ADI( varargin )
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
%argument to this function, the solutio will start from the end of the save
%and rsults will be appended to the previously ran solution.


desired_inputs = 7;

fprintf('Seting up problem to perform Crank-Nicoslon routine.\n');

if nargin < 1, error('Must provide arguments when calling CrankNicolson.m'); end
if nargin ==1,
    if ~isstruct(varargin(1))
        error('Save files must be provided as a structure containing all variables when calling CrankNicolson.n');
    else
        %Load in savefile from previous run
    end
elseif nargin ~= desired_inputs
    error(sprintf('Must provide %i inputs, or one structure containing previous run''s information.',desired_inputs));
else
    %Loading in base variables
    ax = varargin(1);
    bx = varargin(2);
    ay = varargin(3);
    by = varargin(4);
    nodefacx = varargin(5);
    nodefacy = varargin(6);
    DTIMEI = varargin(7);
    
    %Calculate number of internal nodes in X and Y using the following formula
    xnodes = uint32(ceil(exp(ax)));
    ynodes = uint32(ceil(exp(ay)));
    DX = (bx-ax)/(xnodes-2);
    DY = (by-ay)/(ynodes-2);
    x = zeros(xnodes,1);
    y = zeros(ynodes,1);
    for xi=1:xnodes
        x(xi) = DX*xi;
    end
    for yi=1:ynodes
        y(yi) = DY*yi;
    end
    u = zeros(xnodes,ynodes);
    
    %Setting up problem
    TIMEN = 0; %Current time of the problem
    %Setting up guess for initial value of u
    fb = @(x,bx) (bx-x)^2 * cos(pi*x/bx);
    gb = @(x,bx) x*(bx-x)^2;
    for xi=1:xnodes %Top and bottom Dirichlet conditions
        %Assuming u(x,DY) is close to u(x,0)
        u(xi,1) = gb(x(xi),bx);
        %Assuming u(x, by-DY) is close to u(x,by)
        u(xi,ynodes) = fb(x(xi),bx);
        %Initial guess: linear relationship between top and bottom
        deltay = (u(xi,ynodes)-u(xi,1))/(ynodes);
        for yi=2:ynodes-1
            u(xi,yi) = deltay*ui;
        end
        clear deltay
    end
    %Steady state solution time would be quicker if I found an initial
    %condition closer to the steady state solution, but the purpose of this
    %code is to be able to find that steady state solution regardless of
    %the initial condition, so this approximation will do.
    uprev = u; %Solution of previous timestep
end
















