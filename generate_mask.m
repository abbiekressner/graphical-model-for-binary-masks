function mask = generate_mask(ideal,params,verbose,accuracy,max_iter,lookup_table)

  % function mask = generate_mask(ideal,params,verbose,accuracy,max_iter,lookup_table)
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

  %% Organize inputs
  if nargin < 6 || isempty(lookup_table)
    lookup_table = load('example_lookup_table.mat');
  end 
  if nargin < 5
    max_iter = 50;
  end    
  if nargin < 4
    accuracy = 0.01; % (e.g., 0.01 means it'll be accurate to within 1% error)
  end 
  if nargin < 3
    verbose = 1;
  end
  if params.alpha > 1 || params.beta > 1
    params.alpha = params.alpha/100;
    params.beta = params.beta/100;
  end
  limit_params = @(x) max(min(x,0.99),0); % heuristically chosen

  %% Start line search
  [init_A,init_B] = lookup_A_and_B(params);
  params.A = limit_params(init_A);
  params.B = limit_params(init_B);
  iter = 0;
  converged = false;
  while not(converged) && (iter < max_iter)
    iter = iter + 1;
    mask = sample_model(ideal,params,@build_model);
    [~,current_alpha,current_beta] = calc_error(ideal,mask);
    if (abs(current_alpha-params.alpha)<accuracy) && ...
        (abs(current_beta-params.beta)<accuracy)
      converged = true;
    else
      delta = @(x) (0.995/0.25*(-x+0.5).^2)+0.005; % heuristically chosen
      params.A = limit_params( params.A + (params.alpha-current_alpha)*delta(params.A) );
      params.B = limit_params( params.B + (params.beta-current_beta)*delta(params.B) );
    end
    if verbose>0
      fprintf(['Iter %2.0f \t alpha=%3.1f%% (A=%0.4f) \t ' ...
        'beta=%3.1f%% (B=%0.4f) \t gamma=%1.1f\n'], ...
        iter,current_alpha*100,params.A,current_beta*100,params.B,params.gamma);
      if converged
        fprintf('\n');
      end
    end
  end
  if not(converged) && verbose >=0
    fprintf(['\n[WARNING] The mask is inaccurate: gamma= ' num2str(params.gamma,2) ...
      ', alpha= ' num2str(params.alpha,2) ' (actual: ' num2str(current_alpha,2) ')' ...
      ', beta= ' num2str(params.beta,2) ' (actual: ' num2str(current_beta,2) ')\n']);
  end

  function [guess_A,guess_B] = lookup_A_and_B(dummy)

    [val,ind_gamma] = min(abs(lookup_table.g-dummy.gamma));
    grab = @(idx,mat) squeeze(mat(idx,:,:));
    [interp_A,interp_B] = meshgrid(0:0.01:1);
    interp_alphas = interp2(lookup_table.B,lookup_table.A, ...
      grab(ind_gamma,lookup_table.alpha),interp_B,interp_A);
    interp_betas = interp2(lookup_table.B,lookup_table.A, ...
      grab(ind_gamma,lookup_table.beta),interp_B,interp_A);

    [~,ind] = min(abs(interp_alphas(:)-dummy.alpha)+abs(interp_betas(:)-dummy.beta));
    guess_A = interp_A(ind);
    guess_B = interp_B(ind);

  end % end of nested function lookup_A_and_B

end % end of main function
