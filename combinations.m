## Copyright (C) 2023 Leonardo Araujo <leolca@gmail.com>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING. If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{T} = } combinations (@var{A1}, @dots{}, @var{An})
##
## Creates all element combinations of the input arguments @var{A1}, @dots{}, @var{An}.
## Each combination is stored in a row of @var{T}. This function works with numeric 
## and non-numeric arguments but, for non-numeric inputs, it is required to install
## and load the tablicious package.
##
## Example:
## @example
## T = combinations ([1 2 3], [10 20]);
## @end example
##
## @seealso{perms, nchoosek, repelem, repmat}
## @end deftypefn

function T = combinations (varargin)
  tablicious = false;
  use_table = false;
  pkg_list_output = pkg("list");
  for i = 1:length (pkg_list_output), if strcmp (pkg_list_output{i}.name, "tablicious") && pkg_list_output{i}.loaded, tablicious = true; endif endfor
  for i = 1:nargin, is_numeric(i) = isnumeric (varargin{i}); endfor
  all_numeric = all (is_numeric);
  if tablicious,
    use_table = true;
  elseif !all_numeric,
    use_table = true;
    if tablicious,
      warning ("Consider installing tablicious package to get a table output.");
    else,
      error ("Tablicious package is required!");
    endif
  endif

  if use_table,
    variable_names = cell(1, nargin);
    for i = 1:nargin, if !isempty (inputname(i)), variable_names{i} = inputname(i); else variable_names{i} = strcat("Var", num2str (i)); endif; endfor
  endif

  if use_table,
    variables = cell(1, nargin);
    for i = 1:nargin, variables{i} = 1:numel (varargin{i}); endfor
  else
    variables = varargin;
  endif

  for i = 1:nargin, l_args(i) = numel (variables{i}); endfor;
  N = cumprod (fliplr (l_args)); 
  NN = N(end);
  N = fliplr ([1 N(1:end-1)]);
  for i = 1:nargin,
    r = repelem (variables{i}(:), N(i));
    T(:,i) = repmat (r, NN / (l_args(i)*N(i)), 1);
  endfor

  if use_table,
    C = cell(1, nargin);
    for i = 1:nargin, C{i} = varargin{i}(T(:,i)); endfor
    T = cell2table (C, 'VariableNames', variable_names);
  endif
endfunction

%!demo
%! T = combinations ([1 2 3], [10 20])

%!demo
%! T = combinations ([1 2 3], [4 5], [6 7 8 9])

%! ## test input validation
%!error T = combinations ()

%!test
%! A = [1 2 3];
%! B = [4 5];
%! T = combinations (A, B);
%! assert (T, [1 4; 1 5; 2 4; 2 5; 3 4; 3 5])

%!test
%! A = 1;
%! B = [2 3 4 5];
%! C = [6 7];
%! T = combinations (A, B, C);
%! assert (T, [1 2 6; 1 2 7; 1 3 6; 1 3 7; 1 4 6; 1 4 7; 1 5 6; 1 5 7])
