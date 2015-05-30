function [overall,false_positive,false_negative] = calc_error(ideal,mask)

  % function [overall,false_positive,false_negative] = calc_error(ideal,mask)

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

  if all(ismember(unique(mask),[1 2]))
    mask = mask-1;
  end

  overall               = nnz(xor(ideal,mask))/numel(ideal);
  false_positive        = nnz(not(ideal)&mask)/nnz(not(ideal));
  false_negative        = nnz(ideal&not(mask))/nnz(ideal);

  if nargout == 0
    fprintf(['\nOverall error rate: ' repmat(' ',1,5) '%0.1f%%\n'...
      'False positive rate: ' repmat(' ',1,4) '%0.1f%%\n'...
      'False negative rate: ' repmat(' ',1,4) '%0.1f%%\n'],...
      100*overall,100*false_positive,100*false_negative);
  end

end % end of main function
