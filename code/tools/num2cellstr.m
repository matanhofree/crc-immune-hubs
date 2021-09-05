function e=num2cellstr(m)
% numerical matrix to cell str of the same size
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 03, 2011

[r,c]=size(m);
e=cell(r,c);
for i=1:r
    for j=1:c
        e{i,j}=num2str(m(i,j));
    end  
end 
end