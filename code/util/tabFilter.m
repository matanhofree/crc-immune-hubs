function table = tabFilter(x,filterV)
%TABULATE Frequency table.
%   TABLE = TABULATE(X) takes a vector X and returns a matrix, TABLE.
%   The first column of TABLE contains the unique values of X.  The
%   second is the number of instances of each value.  The last column
%   contains the percentage of each value.  If the elements of X are
%   non-negative integers, then the output includes 0 counts for any
%   integers that are between 1 and max(X) but do not appear in X.
%
%   TABLE = TABULATE(X), where X is a categorical variable, character
%   array, or a cell array of strings, returns TABLE as a cell array.  The
%   first column contains the unique string values in X, and the other two
%   columns are as above.
%
%   TABULATE with no output arguments returns a formatted table
%   in the command window.
%
%   See also PARETO.
   
%   Copyright 1993-2011 The MathWorks, Inc.


isnum = isnumeric(x);
if isnum && ~isfloat(x)
    % use of hist() below requires float
    x = double(x);
end
if isnum
   if min(size(x)) > 1,
      error(message('stats:tabulate:InvalidData'));
   end

   y = x(~isnan(x));
else
   y = x;
end

if ~isnum || any(y ~= round(y)) || any(y < 1);
   docell = true;
   [y,yn,yl] = grp2idx(y);
   maxlevels = length(yn);
else
   docell = false;
   maxlevels = max(y);
   %yn = cellstr(num2str((1:maxlevels)'));
end

noSubset = 0;
if nargin<2
    noSubset = 1;
    filterV = true(size(x));
end

%%
[counts, values] = hist(y,(1:maxlevels));
[countsF, valuesF] = hist(y(filterV),(1:maxlevels));
%%
total = sum(counts);
%percents = 100*counts./total;
percents = 100*countsF./counts;

if exist('yn','var')
    [~,zidxOrd] = sort(yn);
else
    zidxOrd = 1:maxlevels;
end

%%

if nargout == 0
   if docell
      width = max(cellfun('length',yn));
      width = max(5, min(50, width));
   else
      width = 5;
   end
   
   if noSubset
       % Create format strings similar to:   '  %5s    %5d    %6.2f%%\n'
       fmt1 = sprintf('  %%%ds    %%6d    %%6.2f%%%%\n',width);
       fmt2 = sprintf('  %%%ds    %%6s    %%6s\n',width);
       fmt3 = sprintf('--%s\n  %%%ds    %%6d    %%6.2f%%%%\n',repmat('-',1,width'),width);
       fprintf(1,fmt2,'Value','Count','Percent');
       if docell
          for zi=1:maxlevels
             j = zidxOrd(zi);
             fprintf(1,fmt1,yn{j},countsF(j),100*countsF(j)/total);         
          end
          fprintf(1,fmt3,'Total:',sum(countsF),100*sum(countsF)/total);  
       else
          fprintf(1,'  %5d    %5d    %6.2f%%\n',[values' counts' percents']');
       end              
   else
       % Create format strings similar to:   '  %5s    %5d    %6.2f%%\n'
       fmt1 = sprintf('  %%%ds    %%6d    %%5d    %%6.2f%%%%\n',width);
       fmt2 = sprintf('  %%%ds    %%6s    %%5s    %%6s\n',width);
       fmt3 = sprintf('--%s\n  %%%ds    %%6d    %%5d    %%6.2f%%%%\n',repmat('-',1,width'),width);
       fprintf(1,fmt2,'Value','Filter','Orig.','Percent');
       if docell
          for zi=1:maxlevels
             j = zidxOrd(zi);
             fprintf(1,fmt1,yn{j},countsF(j),counts(j),percents(j));         
          end
          fprintf(1,fmt3,'Total:',sum(countsF),total,100*sum(countsF)/total);  
       else
          fprintf(1,'  %5d    %5d    %6.2f%%\n',[values' counts' percents']');
       end
   end
else
   if ~docell
      table = [values' countsF' counts' percents'];
   elseif isnum
      table = [yl(:) countsF' counts' percents'];
   else
      table = [yn num2cell([countsF' counts' percents'])];
   end
   
   table = table(zidxOrd,:);
end
