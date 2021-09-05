function [b, m, n, cnt,dup_list] = uniquec(A)
% A unique function that doesn't break the ordering things are given in
% Outputs:
% b unique list
% b = A(m)
% A = b(n)
% cnt = count(b(i)))
% dup_list = list of occurances of each item in original array!



    m = [];
    n = [];
    
    z = 1;
    
    if (iscell(A))
        pos_vector = containers.Map('keyType','char','valueType','double');
        pos_vector(A{1}) = z;
        dup_list{z} = z;
        cnt(z) = 1;
        
        m(z) = 1;
        n(1) = 1;
        
        ct = 2;
        for i = 2:length(A)
            if ( pos_vector.isKey(A(i)))
                z = pos_vector(A{i});
                dup_list{z} = [ dup_list{z} i ];
                cnt(z) = cnt(z) + 1;
                n(i) = z;
            else
                z = ct;
                pos_vector(A{i}) = z;
                dup_list{z} = i;
                cnt(z) = 1;
                m(z) = i;
                n(i) = z;
                
                ct = ct + 1;
            end
        end
    else
        pos_vector = containers.Map('keyType','double','valueType','double');
        pos_vector(A(1)) = z;
        dup_list{z} = z;
        cnt(z) = 1;
        m(z) = 1;
        n(1) = 1;
          
        ct = 2;
        for i = 2:length(A)
            if ( pos_vector.isKey(A(i)))
                z = pos_vector(A(i));
                cnt(z) = cnt(z) + 1;
                dup_list{z}(cnt(z)) = i;                
                n(i) = z;
            else
                z = ct;
                pos_vector(A(i)) = z;
                dup_list{z} = i;
                cnt(z) = 1;
                m(z) = i;
                n(i) = z;
                
                ct = ct + 1;                                
            end
        end        
    end
    
    [~,oidx] = sort(cell2mat(pos_vector.values()));
    b = pos_vector.keys();
    b = b(oidx);
    
    if (~iscell(A))
        b= cell2mat(b);
    end

end