
function [elapsed_time] = funNS_solver_dt(const_dt, N)

    % a function running the iterations returning the time it took to converge
    
    warning off
    clearvars -except const_dt N


    Re = 1000;             % Reynolds number
%     N = 48;                % Number of volumes in the x- and y-direction
    Delta = 1/N;           % uniforme ruimte stap in de afbeelding

%     dt = const_dt/(N^4);
    dt = const_dt;
    % dt = 0.000001;           % time step




    tol = 1e-6;            % tol determines when steady state is reached and the program terminates

    % wall velocities
    U_wall_top = -1;
    U_wall_bot = 0;
    U_wall_left = 0;
    U_wall_right = 0;
    V_wall_top = 0;
    V_wall_bot = 0;
    V_wall_left = 0;
    V_wall_right = 0;

    %
    %   Generation of a non-uniform mesh
    %

    %
    %   x are the coordinates of the nodal points on the outer-oriented grid
    %
    tx = zeros(1,N+1);
    for i=1:N+1
        xi = (i-1)*Delta;
        tx(i) = 0.5*(1. - cos(pi*xi));
    end

    % Local mesh size on outer oriented grid
    th = tx(2:N+1) - tx(1:N);


    %
    %  x are the coordinates of the nodal points on the inner-oriented grid (including
    %  endpoints 0 and 1)
    %  h contains the edge lengths on the inner-oriented grid
    %
    x = 0.5*(tx(1:N) + tx(2:N+1));
    x = [0 x 1];

    h = x(2:N+2) - x(1:N+1);

    %
    %   Initial condition u=v=0
    %
    %   Both u and v will be stored in one big vector called 'u'
    %
    %   The vector u only contains the true unknowns, not the velocities that
    %   are prescribed by the boundary conditions
    %
    %   The vector u contains the *inner-oriented* fluxes as unknowns
    %
    u = zeros(2*N*(N-1),1);

    % Set up the Incindence matrix 'tE21' which connects the fluxes to the
    % volumes. Use the orientation described in the assignment.

    ............... % DONE
    s = repmat([-1 1 -1 1], 1, N*N);
    i = zeros(size(s));
    j = zeros(size(s));
    for ii = 1:N*N
        i((1:4) + (ii-1)*4) = [ii ii ii ii];
        j((1:4) + (ii-1)*4) = [ii+floor((ii-1)/N)  ii+1+floor((ii-1)/N)  N*(N+1)+ii  N*(N+1)+N+ii];
    end

    tE21 =sparse(i,j,s);

    % full(tE21)

    %
    %  Inserting boundary conditions for normal velocity components
    % 
    ................ % DONE
    % first find the columns of tE21 which contain only one value as these are the boundaries
    [~, columns_tE21, values_tE21] = find(tE21); 
    col_diff_tE21 = diff(columns_tE21);
    col_diff_tE21(circshift(col_diff_tE21, 1)==0) = 0;
    col_diff_tE21(end+1) = col_diff_tE21(end);                                  % col_diff contains a 1 where the collumns are with only one value
    col_idx_tE21 = columns_tE21(col_diff_tE21 == 1);                            % indices of the columns containing only one value

    u_norm_e = tE21(:, col_idx_tE21);                                % orientation of the normal velocities
    u_norm_bc = [repmat([U_wall_left U_wall_right], 1, N)...
                V_wall_bot*ones(1,N) V_wall_top*ones(1,N)];                     % magnitude of normal velocities at the wall
    u_norm = u_norm_e*u_norm_bc';                                               % multiply direction and magnitude to get u norm vector 

    %
    %   Remove columns associated with prescribed normal velocities from
    %   Incidence matrix 'tE21'
    %
    ................. % DONE
    tE21(:, col_idx_tE21) = [];                                                 % remove columns with those indices


    % Setting up simple Hodge matrix which converts the fluxes on the
    % outer-oriented grid to circulation on the inner-oriented grid. Assume
    % that the fluxes and circulation are constant over each 1-cell. This will
    % give a diagonal Hodge matrix. Call this Hodge matrix 'H1t1'
    .................% DONE
    tsigma = zeros(2*N*(N+1),1);
    sigma = zeros(2*N*(N+1),1);

    % make tsigma vector containing the primal grid lengths where the fluxes act on
    for ii = 1:N
        tsigma((1:N+1)+(ii-1)*(N+1)) = th(ii);
    end
    tsigma(1+N*(N+1):end) = repmat(th, 1 ,N+1);

    % make sigma vector containing dual grid lengths where the velocities act along
    sigma(1:N*(N+1)) = repmat(h, 1, N);
    for ii = 1:N+1
        sigma((1:N)+(ii-1)*(N)+ N*(N+1)) = h(ii);
    end

    % divide those to find diagonal of hodge matrix, create hodge matrix
    H1t1 = sparse(1:length(sigma),1:length(sigma),sigma./tsigma);


    % Hu_norm = zeros(2*N*(N+1),1); Dont need this

    %
    % Hu_norm is the vector which will contain the Hodge of the prescribed
    % normal fluxes. Calculate the vector 'Hu_norm' (H*u_norm I presume)
    ............... % DONE
    H1t1_norm = H1t1(col_idx_tE21,:);
    H1t1_norm = H1t1_norm(:,col_idx_tE21); % remove empty rows
    % Hu_norm = H1t1_norm*u_norm'; 

    %
    %  Remove corresponding row and columns from the Hodge matrix and also
    %  Remove the corresponding 'rows' from Hu_norm (already done)
    ................. % DONE
    H1t1(col_idx_tE21,:) = [];
    H1t1(:,col_idx_tE21) = [];

    %
    % Set up the incidence E^{2,1} between 1-cochain circulation and 2-cochain vorticity on
    % the inner-oriented (extended) grid
    %
    % This incidence matrix will be called 'E21' in the program
    .................... % DONE
    s = repmat([1 -1 -1 1], 1, (N+1)*(N+1));
    i = zeros(size(s));
    j = zeros(size(s));

    for ii = 1:(N+1)*(N+1)
        i((1:4) + (ii-1)*4) = [ii ii ii ii];
        j((1:4) + (ii-1)*4) = [ii  ii+N+1 (N+2)*(N+1)+ii+floor((ii-1)/(N+1))  (N+2)*(N+1)+ii+floor((ii-1)/(N+1))+1];
    end

    E21 =sparse(i,j,s);

    % Inserting prescribed tangential bundary conditions
    ....................% DONE

    % first find the columns of E21 which contain only one value as these are the boundaries
    [~, columns_E21, values_E21] = find(E21); 
    col_diff_E21 = diff(columns_E21);
    col_diff_E21(circshift(col_diff_E21, 1)==0) = 0;
    col_diff_E21(end+1) = 1;                                                    % add a 1 on the end because the last value is a boundary
    col_idx_E21 = columns_E21(col_diff_E21 == 1);                               % indices of the columns containing only one value in E21 


    u_tang_bc = [U_wall_bot*h U_wall_top*h...
                 repmat([V_wall_left V_wall_right], 1, N+1)];                   % magnitude of tangential velocities at the wall

    % Find adjusted column indices of the known normal velocities
    col_idx_tE21_adj = (1:max(size(E21)))';
    col_idx_tE21_adj(col_idx_E21) = [];
    col_idx_tE21_adj = col_idx_tE21_adj(col_idx_tE21);

    % Join normal and tangential column indices and boundary conditions
    u_pres_idx = [col_idx_tE21_adj; col_idx_E21];
    u_pres_bc = [u_norm_bc'; u_tang_bc'];

    % Sort joined column indices and boundary conditions
    [u_pres_idx, sorted_idx] = sort(u_pres_idx); 
    u_pres_bc = u_pres_bc(sorted_idx); % contains the magnitude of the boundary conditions for the known normal and tangential velocities

    u_pres_e = E21(:,u_pres_idx); % contains the indicence matrix for the known normal and tangential velocities

    % Remove columns from the incidence matrix E21 corresponding to both the
    % prescribed tangental velocities and normal velocities
    ................. % DONE
    E21(:,u_pres_idx) = []; 

    %
    % Store the prescribed normal and tangential velcoities in the vector
    % 'u_pres'
    ........................ % DONE
    % u_pres = u_pres_e*u_pres_bc; % ?? change here
    u_pres = u_pres_bc;

    % Set up the Hodge matrix which maps inner-oriented 2-cochains to
    % outer-oriented 0-cochains. Assume that the vorticity is constant in the
    % inner-oriented 2-cells. This will give a diagonal Hodge matrix. Call this
    % Hodge matrix 'Ht02'
    ...................
    Ht02_vec = zeros((N+1)*(N+1), 1);
    for ii = 1:length(h)
        Ht02_vec((1:(N+1))+ (ii-1)*(N+1)) = 1./(h(ii)*h);
    end

    Ht02 = sparse(1:length(Ht02_vec),1:length(Ht02_vec),Ht02_vec);


    % Set up the Hodge matrix which maps inner-oriented 1-cochains to
    % outer-oriented 1-cochains. Call this Hodge matrix 'Ht11'. Assume again
    % that everything is constant over the 1-cells, which will then yield a
    % diagonal Hodge matrix.
    ..................
    % Ht11 is the inverse of H1t1
    Ht11 = H1t1^-1;


    % The prescribed velocties will play a role in the momentum equation


    % u_pres = H1t1*E21'*Ht02*u_pres; ?? does this need to go? is in RHS
    % poisson eq 


    % 
    % Now all matrices are set up and the time stepping cam start. 'iter' will
    % record the number of time steps. This allows you to give output after a
    % preselected number of time steps.
    % 
    % 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
    % small, steady state has been reached. Determine a suitable value for
    % 'tol'


    difff = 1;
    iter = 1;
    diff_list = [];

    % Create i and j lists for the u and v convective terms: (same for each iteration)

    i_u = repmat(2:N, 1, N);
    j_u = ceil((1:N*(N-1))/(N-1));

    i_v = repmat(1:N, 1, N-1);
    j_v = ceil((N+1:N+N*(N-1))/(N));

    E21_pres = u_pres_e;

    % Set up the matrix for the POisson equation
    A = -tE21*Ht11*tE21'; % Corrected H1t1 to Ht11 
    tic

    while difff > tol
    %    disp(iter)


        %
        %   Calculate the convective terms using the vorticty and the local
        %   velocity field. Store the convective terms in a vector
        %   'convective'.
        %
        %  Note that you only have to set up the convective terms for true
        %  velicoty unknowns and not for those inner-oriented circulations for
        %  which the value is already known.
        %   
        ...................... % DONE

        xi = E21*u + E21_pres*u_pres;   % xi on inner oriented grid with U-pres
        txi = Ht02*xi;                  % txi on outer oriented grid with U-pres

        % add prescribed values to u vector to speed up calculation of convective terms
        idx_list = 1:2*N*(N+1);
        u_full = zeros(2*N*(N+1),1);
        u_full(setdiff(idx_list, col_idx_tE21)) = u;


        conv_u = -h(i_u)./(4*h(j_u)).*(findv(i_u-1,j_u, N, u_full, col_idx_tE21) + findv(i_u,j_u, N, u_full, col_idx_tE21)).*findxi(i_u,j_u, N, txi)...
                -h(i_u)./(4*h(j_u+1)).*(findv(i_u-1,j_u+1, N, u_full, col_idx_tE21) + findv(i_u,j_u+1, N, u_full, col_idx_tE21)).*findxi(i_u,j_u+1, N, txi);

        conv_v = h(j_v)./(4*h(i_v)).*(findu(i_v,j_v-1, N, u_full, col_idx_tE21) + findu(i_v,j_v, N, u_full, col_idx_tE21)).*findxi(i_v,j_v, N, txi)...
                +h(j_v)./(4*h(i_v+1)).*(findu(i_v+1,j_v-1, N, u_full, col_idx_tE21) + findu(i_v+1,j_v, N, u_full, col_idx_tE21)).*findxi(i_v+1,j_v, N, txi);

        convective = [conv_u conv_v]';


        % Set up the right hand side for the Poisson equation for the pressure


        rhs_Poisson  =   tE21*Ht11*(u/dt  - convective - H1t1*E21'*Ht02*(xi/Re))+ u_norm/dt; 
        % xi is calculated before so the calculation here is removed 
        % changed Ht11*E21'*Ht02* to H1t1*E21'*Ht02* as this was wrong 

        % calculation for A moved outside iteration as this does not change

        % Solve for the new pressure

        P = A\rhs_Poisson;

        % Store the velocity from the previous time step in the vector u_old

        uold = u;

        % Udate the velocity field

        u = u - dt* (convective - tE21'*P + H1t1*E21'*Ht02*(xi/Re));

        %
        %  Every other 1000 iterations check whether you approach steady state
        %  and check whether you satisfy conservation of mass. The largest
        %  rate at which mass is destroyed or created is denoted by 'maxdiv'.
        %  This number should be very small, in the order of machine precision.


        maxdiv = max(tE21*Ht11*u + u_norm);

        difff = max(abs(u-uold))/dt;

        
%         diff = 1e-7 % used in debugging 
        if mod(iter,1000) == 0
            disp(difff)
        end  

%         figure(1)
%         plot(iter/1000, difff, 'bx')
%         hold on

        diff_list(end+1) = difff;
        elapsed_time = toc;
        

%         if elapsed_time > 300
%             break
%         end
        if isnan(difff)
            elapsed_time = 10000;
        else
        if difff > 0.5 && iter > 1000
            disp('diff above 0.5 after 500 iter')
            elapsed_time = 10000;
            break
        end
            
        iter = iter + 1;
    end
%     disp(elapsed_time)
end






