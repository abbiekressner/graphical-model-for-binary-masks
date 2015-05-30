function masks = sample_model(ideal_mask,params,build_func,nSamples,burnIn,init);

  % function masks = sample_model(ideal_mask,params,build_func,nSamples,burnIn,init);
  %
  % ideal_mask: ideal binary mask to build from
  %
  % params: structure containing parameters for the build_func
  %
  % build_func: handle to function that builds the model [default @build_model]
  %
  % nSamples: number of masks to get [default 1]
  %
  % burnIn:   number of Gibbs "burn-in" samples [default 100
  % (heuristically chosen)]
  %
  % init:     initial states of unobserved nodes [defaults to IBM]
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
  
  %% Initialize variables
  observed = ideal_mask + 1; % make 0s and 1s into indices of states
  [nChans,nFrames] = size(observed);
  nTimeFreqUnits = nChans*nFrames;

  %% Set defaults
  if nargin < 6
    init = reshape(observed,[],1); % init as IBM (or max likelihood?)
  end
  if nargin < 5
    burnIn = 100;
  end
  if nargin < 4
    nSamples = 1;
  end
  if nargin < 3
    build_func = @build_model;
  end

  %% Build the UGM
  [nodePot,edgePot,edgeStruct,clamped] = build_func(ideal_mask,params);

  %% Get sample(s)
  edgeStruct.maxIter = nSamples;
  samples = UGM_Sample_Conditional(nodePot,edgePot,edgeStruct,clamped,...
    @(x,y,z) UGM_Sample_Gibbs(x,y,z,burnIn,init));

  %% Reshape the bm matrix
  masks = NaN(nChans,nFrames,nSamples);
  for ii = 1:nSamples
    masks(:,:,ii) = reshape(samples(1:nTimeFreqUnits,ii),...
      [nChans nFrames]);
  end
  masks = masks - 1; % make indices of states into 0s and 1s again

end % end of main function
