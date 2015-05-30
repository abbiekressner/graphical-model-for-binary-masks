function lookup_table = build_lookup_table(eg_ideal,g,verbose)

  if nargin < 3
    verbose = true;
  end
  if nargin < 2
    g = 1.0:0.1:3.5;
  end
  A = [0:0.05:0.4    0.425:0.025:0.60    0.7:0.1:0.9    0.99];
  B = A;
  
  alpha = NaN(numel(g),numel(A),numel(B));
  beta = NaN(numel(g),numel(A),numel(B));

  cnt = 0;
  if verbose, textprogressbar('Building the lookup table... '); end;
  for xx = 1:numel(g)
    for yy = 1:numel(A)
      for zz = 1:numel(B)
        [~,alpha(xx,yy,zz),beta(xx,yy,zz)] = ...
          calc_error(eg_ideal,sample_model(eg_ideal, ...
          struct('gamma',g(xx),'A',A(yy),'B',B(zz))));
        cnt = cnt + 1;
        if verbose, textprogressbar(100*cnt/numel(alpha)); end;
      end
    end
  end
  if verbose, textprogressbar('Done.'); end;

  lookup_table.g = g;
  lookup_table.A = A;
  lookup_table.B = B;
  lookup_table.alpha = alpha;
  lookup_table.beta = beta;

end % end of main function
