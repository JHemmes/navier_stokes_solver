function [out] = findxi(i_in,j_in, N, xi)

    % returns Xi_{i,j} for i_in, j_in. i_in and j_in can be lists

    out = zeros(size(i_in));

    for kk = 1:length(i_in)
        i = i_in(kk);
        j = j_in(kk);
        index = i+(j-1)*(N+1);
        out(kk) = xi(index);
    end
    
end
