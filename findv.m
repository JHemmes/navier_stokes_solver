function [out] = findv(i_in,j_in, N, u, pres_col)

    % returns v_{i,j} for i_in, j_in. i_in and j_in can be lists
    % u can be the vector with true unknowns, or vector with also the known prescribed values. 
    % Adding the known values to u omits the timeconsuming ismember operation

    out = zeros(size(i_in));
    
    if length(u) == 2*N*(N+1)
        for kk = 1:length(i_in)
            i = i_in(kk);
            j = j_in(kk);
            index = N*(N+1) + i + (j-1)*N;
            out(kk) = u(index);  
        end

        
    else
        for kk = 1:length(i_in)
            i = i_in(kk);
            j = j_in(kk);
            index = N*(N+1) + i + (j-1)*N;

            if ismember(index, pres_col)
                out(kk) = 0; % note that this should actually come from u_norm_bc. This way the normal velocity BCs are hardcoded to be 0
            else 
                idx_list = 1:2*N*(N+1);
                idx_list(pres_col) = [];
                out(kk) = u(idx_list == index);
            end
        end
    end
    
end