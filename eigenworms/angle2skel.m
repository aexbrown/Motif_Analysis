function [skelX, skelY] = angle2skel(angleArray, meanAngle, arclength)

% ANGLE2SKEL Take in an angle array and integrate over the angles to get
% back a skeleton for each frame.  NB: This reconstruction assumes each 
% segment was equally spaced so that each reconstructed skeleton segment 
% has length arclength/(numAngles + 1)
% 
% Input
%   angleArray - a numSkelPoints - 1 by numFrames array of skeleton
%                tangent angles that have been rotated to have a mean
%                angle of zero.
%   meanAngle  - a 1 by numFrames array of angles.  Each angle is the mean
%                angle of the skeleton used to make the corresponding row
%                of angles in angleArray.  Can be left as zeros if no
%                rotation is desired.
%   arclength  - the total arclength of the skeleton to be reconstructed.
%                Can be set to 1 for a normalised skeleton.
% 
% Output
%   skelX      - a numAngles + 1 by numFrames array of skeleton
%                x-coordinates
%   skelY      - a numAngles + 1 by numFrames array of skeleton
%                y-coordinates
% 
% 
% Copyright Medical Research Council 2013
% André Brown, andre.brown@csc.mrc.ac.uk, aexbrown@gmail.com
% 
% 
% The MIT License
% 
% Copyright (c)  Medical Research Council 2013
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


% get dimensions
numAngles = size(angleArray, 1);
numFrames = size(angleArray, 2);


% initialisation
skelX = NaN(numAngles + 1, numFrames);
skelY = NaN(numAngles + 1, numFrames);

for ii = 1:numFrames
    % add up x-contributions of angleArray, rotated by meanAngle
    skelX(:, ii) = [0; cumsum(cos( angleArray(:, ii) + meanAngle(ii) ) * ...
        arclength/(numAngles + 1)) ];
    
    % add up y-contributions of angleArray, rotated by meanAngle
    skelY(:, ii) = [0; cumsum(sin( angleArray(:, ii) + meanAngle(ii) ) * ...
        arclength/(numAngles + 1)) ];
end


