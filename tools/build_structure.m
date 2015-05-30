function [edgeStruct,nNodes] = build_structure(nChans,nFrames,nLayers,useMex)

  % TODO fill in the help info
  %
  % --------------------------------------------------------------------
  % Written by Abigail Kressner in 2014.
  %
  % This file is part of GMBM.
  %
  % GMBM is free software: you can redistribute it and/or modify
  % it under the terms of the GNU General Public License as published by
  % the Free Software Foundation, either version 3 of the License, or
  % (at your option) any later version.
  %
  % GMBM is distributed in the hope that it will be useful,
  % but WITHOUT ANY WARRANTY; without even the implied warranty of
  % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  % GNU General Public License for more details.
  %
  % You should have received a copy of the GNU General Public License
  % along with GMBM.  If not, see <http://www.gnu.org/licenses/>.
  % --------------------------------------------------------------------

% Initialize
if nargin < 3
  nLayers = 1;
end
if nargin < 4
  useMex = true;
end
nTimeFreqUnits = nChans*nFrames;
nNodes = nLayers*nTimeFreqUnits;
adj = sparse(nNodes,nNodes);

% Add Down Edges
ind = 1:nTimeFreqUnits;
exclude = sub2ind([nChans nFrames],...
  repmat(nChans,[1 nFrames]),1:nFrames); % no Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

% Add Right Edges
ind = 1:nTimeFreqUnits;
exclude = sub2ind([nChans nFrames],...
  1:nChans,repmat(nFrames,[1 nChans])); % no right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nChans)) = 1;

% Add Up/Left Edges
adj = adj+adj';

if nLayers > 1

  % Connect base layer to top ("observed") layer
  ind = 1:nTimeFreqUnits;
  adj(sub2ind([nNodes nNodes],ind,ind+nTimeFreqUnits)) = 1;
  adj(sub2ind([nNodes nNodes],ind+nTimeFreqUnits,ind)) = 1;

end

% Call makeEgdeStruct
nStates = 2; % binary mask has only two states
edgeStruct = UGM_makeEdgeStruct(adj,nStates,useMex);
