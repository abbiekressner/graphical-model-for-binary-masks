function [gamma,model,potentials_and_marginals] = learn_parameters(ideal,nonideal,verbose,optimization)

  % function [gamma,model,potentials_and_marginals] = ...
  %   learn_parameters(ideal,nonideal,verbose,optimization)
  %
  % ideal: cell array containing ideal binary masks
  % nonideal: cell array containing corresponding non-ideal binary masks
  %
  % model: struct containing the model
  % potentials_and_marginals: struct containing node and edge potentials
  %
  % Note: gamma (as defined in Kressner and Rozell, JASA, 2015) can be seen
  % along the diagonals in potentials_and_marginals.edgePot(:,:,1). Also,
  % by definition potentials_and_marginals.edgePot(:,:,1) is the same as
  % potentials_and_marginals.edgePot(:,:,2) and as
  % potentials_and_marginals.edgePot(:,:,3) and so on.
  %
  % See http://www.cs.ubc.ca/~schmidtm/Software/UGM/trainCRF.html and
  % http://www.cs.ubc.ca/~schmidtm/Software/UGM/trainMRF.html for guidance
  % through the following code. In particular, see the "Using Node Features"
  % section of the trainCRF page.
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

if nargin < 4
  optimization = 'l2';
end
if nargin < 3
  verbose = true;
end


%% Initialize and truncate masks
useMex = true;
nChans = size(ideal{1},1);
minFrames = min(cellfun(@(x) size(x,2),ideal)); % TODO incorporate UGM_CRFcell_NLL so that I don't have to truncate the masks 
nInstances = numel(ideal);
nNodes = nChans*minFrames; % = nTimeFreqUnits
features = zeros(nInstances,nNodes,'int32');
labels = zeros(nInstances,nNodes,'int32');
for ii = 1 : nInstances 
  % Make 0 and 1 into indices of states (1 and 2)
  features(ii,:) = reshape(...
    ideal{ii}(:,1:minFrames)+1,   1,[]);
  labels(ii,:) = reshape(...
    nonideal{ii}(:,1:minFrames)+1,   1,[]);  
end
nStates = 2;
if verbose, fprintf('\n\nInitialization complete.\n'); end;


%% Build edge structure
nLayers = 1;
edgeStruct = build_structure(nChans,minFrames,nLayers,useMex);
if verbose, fprintf('Built edgeStruct.\n'); end;


%% Build Xnode
%Xnode(1:nInstances,1,1:nNodes) = double(features == 2);
nFeatures = 2; % binary indicator variable for what dominates 
Xnode = zeros(nInstances,nFeatures,nNodes);
for ff = 1 : nFeatures
  [idx1,idx2] = find(features==ff);
  for ii = 1 : numel(idx1)
    Xnode(idx1(ii),ff,idx2(ii)) = 1;
  end
end
Xnode = [ones(nInstances,1,nNodes) Xnode]; % add bias variable
nNodeFeatures = size(Xnode,2);
if verbose, fprintf('Built Xnode.\n'); end;


%% Build nodeMap
nodeMap = zeros(nNodes,nStates,nNodeFeatures,'int32');
nodeMap(:,2,1) = 1; % bias
nodeMap(:,1,2) = 2;
nodeMap(:,2,3) = 3;
if verbose, fprintf('Built nodeMap.\n'); end;


%% Build Xedge (note: edges are only among unobserved nodes here)
% Option 1 // bias variable (not dependent on ideal mask)
nEdges = edgeStruct.nEdges;
Xedge = ones(nInstances,1,nEdges);

% Option 2 // dependent on nodes of the ideal mask (aka features for
% edges)
%sharedFeatures = [1; zeros(nFeatures,1)];
%Xedge = UGM_makeEdgeFeatures(Xnode,edgeStruct.edgeEnds,sharedFeatures);

if verbose, fprintf('Built Xedge.\n'); end;


%% Build edgeMap
% Option 1 // same state vs different state
nEdges = edgeStruct.nEdges;
edgeMap = zeros(nStates,nStates,nEdges,'int32');
edgeMap(1,1,:) = max(nodeMap(:)) + 1;
edgeMap(2,2,:) = max(nodeMap(:)) + 1;
% That is, we will make all edges use the n-th element of the parameter
% vector w for the potential of adjacent TF units being in the same
% state. And further, we will fix the potential of adjacent TF units
% having different states.

% Option 2 // different params for all different combinations of states
%nEdges = edgeStruct.nEdges;
%edgeMap = zeros(nStates,nStates,nEdges,'int32');
%edgeMap(1,1,:) = max(nodeMap(:))+1;
%edgeMap(2,1,:) = max(nodeMap(:))+2;
%edgeMap(1,2,:) = max(nodeMap(:))+3;

% Option 3 // map "edge features" (corresponds to Option 2 for Xedge)
%nEdges = edgeStruct.nEdges;
%nEdgeFeatures = size(Xedge,2);
%edgeMap = zeros(nStates,nStates,nEdges,nEdgeFeatures,'int32');
%ff = max(nodeMap(:));
%for edgeFeat = 1:nEdgeFeatures
  %for s1 = 1:nStates
    %for s2 = 1:nStates
      %ff = ff + 1;
      %edgeMap(s1,s2,:,edgeFeat) = ff;
    %end
  %end
%end

if verbose, fprintf('Built edgeMap.\n'); end;


%% Learn parameters
if verbose, fprintf('Starting to learn the parameters...\n\n'); end;

if strcmp(optimization,'ml')
  % Option 1 // maximum likelihood estimation
  nParams = max([nodeMap(:);edgeMap(:)]);
  w = zeros(nParams,1);
  if verbose
    options = [];
  else
    options = struct('Display','off');
  end
  w = minFunc(@UGM_CRF_NLL,w,options,...
    Xnode,Xedge,labels,nodeMap,edgeMap,edgeStruct,...
    @UGM_Infer_MeanField);
elseif strcmp(optimization,'l2')
  % Options 2 // L2-regularized parameter estimation
  nParams = max([nodeMap(:);edgeMap(:)]);
  lambda = 10*ones(nParams,1); % TODO which lambda values are best?
  lambda(1) = 0; % don't penalize bias
  regFunObj = @(w) penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,... 
    labels,nodeMap,edgeMap,edgeStruct,@UGM_Infer_MeanField);
  w = zeros(nParams,1);
  if verbose
    options = [];
  else
    options = struct('Display','off');
  end
  w = minFunc(regFunObj,w,options);
else
  error('Armageddon!!!');
end


%% Gather outputs
model.w = w;
model.Xnode = Xnode;
model.Xedge = Xedge;
model.nodeMap = nodeMap;
model.edgeMap = edgeMap;
model.edgeStruct = edgeStruct;

%% Make results readable
if nargout > 0
  [potentials_and_marginals.nodePot,potentials_and_marginals.edgePot] = ...
    UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct);
  [potentials_and_marginals.nodeMar,potentials_and_marginals.edgeMar] = ...
    UGM_Infer_MeanField(potentials_and_marginals.nodePot,...
    potentials_and_marginals.edgePot,edgeStruct);
  gamma = potentials_and_marginals.edgePot(1,1,1);  
end
