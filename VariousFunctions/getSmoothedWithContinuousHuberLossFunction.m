%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth an if-statement by using a Huber loss function. 
%
% if-statements in the calculation of the Huber loss function are
% approximated with tanh, to make the result compatible with algorithmic
% differentiation. For sufficiently high scale factors inside the
% tanh-function, this approximation is exact up to machine precision.
%
% Adapted from: https://github.com/opensim-org/opensim-core/blob/master/
%       OpenSim/Simulation/Model/Bhargava2004SmoothedMuscleMetabolics.cpp#L188
%
% Input arguments
%       cond: condition
%       left: ideal output if condition < 0
%       right: ideal output if condition > 0
%       smoothing: factor to determine how smooth the transition between
%           left and right is. Lower value is smoother, higher value is
%           more accurate. (passed in varargin{1})
%       direction: unsure, but having it = -1 works (passed in varargin{2})
%
% Output argument
%       ret: return value
%
% Author: Lars D'Hondt
% Date: 20/Dec/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret = getSmoothedWithContinuousHuberLossFunction(cond,left,right,varargin)

if length(varargin)>=1
    smoothing = varargin{1};
else
    smoothing = 10;
end
if length(varargin)>=2
    direction = varargin{2};
else
    direction = -1;
end

offset = left*ifPos(direction) + right*ifNeg(direction);

scale = (right-left)/cond;
delta = 1;
state = direction*cond;
shift = 0.5*(1/smoothing);
y = smoothing*(state+shift);

f = offset + (0.5.*y.*y).*ifPos(y).*ifNeg(y-delta) + (delta.*(y-0.5*delta)).*ifPos(y-delta);

ret = scale*(f/smoothing + offset*(1-1/smoothing));

function factor = ifPos(cond)
    factor = 0.5 + 0.5*tanh(cond*1e6);
end
function factor = ifNeg(cond)
    factor = 0.5 - 0.5*tanh(cond*1e6);
end
end