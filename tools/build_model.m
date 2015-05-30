function [nodePot,edgePot,edgeStruct,clamped] = build_model(ideal,params)

  % function [nodePot,edgePot,edgeStruct,clamped] = build_model(ideal,params)
  %
  % ideal: ideal binary mask to build from  
  %
  % params: structure containing the following fields
  %
  % gamma: neighbors in the masks are this much more
  % likely to take the same state than they are to take different states
  %
  % A: likelihood factor that the masker-dominated units are
  % labeled target-dominated
  %
  % B: likelihood factor that the target-dominated units are
  % labeled masker-dominated 
  %
  % lambda: prior probability of being in each state
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

  if not(isfield(params,'lambda')) 
    params.lambda = 0.5; % potential of being masker-dominated
  end    

  %% Set parameters 
  observed = ideal + 1; % make 0s and 1s into indices of states
  useMex = true;

  %% Initialize
  [nChans,nFrames] = size(observed);
  nTimeFreqUnits = nChans*nFrames;
  nLayers = 2; % observed and unobserved
  nStates = 2;

  %% Build edge structure
  [edgeStruct,nNodes] = build_structure(nChans,nFrames,nLayers,useMex);

  %% Build node potentials
  nodePot = zeros(nNodes,nStates);
  nodePot(:,1) = params.lambda;
  nodePot(:,2) = 1-params.lambda;

  %% Set edge potentials 
  edgePot = zeros(nStates,nStates,edgeStruct.nEdges);
  for ee = 1:edgeStruct.nEdges
    if all(edgeStruct.edgeEnds(ee,:) <= nTimeFreqUnits)
      % edge is between neighboring unobserved nodes
      edgePot(:,:,ee) = [params.gamma 1;
                        1 params.gamma];
    else 
      % edge is between unobserved and observed
      edgePot(:,:,ee) = [1-params.A params.B;
                        params.A 1-params.B]; 
    end
  end

  %% Declare observations (via "clamping")
  clamped = zeros(nNodes,1);
  clamped(nTimeFreqUnits+1:nNodes) = observed(:);

end % end of main function
