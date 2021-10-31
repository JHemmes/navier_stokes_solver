function [totaldiff, elapsed_time, iter] = funNS_solver_tol(tol, N)
    
    % function used to determine the error, elapsed time and number of
    % iterations used to find the optimal stoppng tolerance
    

    warning off

    % This file contains the skeleton of the program which will solve the lid
    % driven cavity problem on a unit square. The parts that have to be
    % supplemented are described in the assignment.
    %
    % The pieces that need to be filled in are indicated by dots: ............
    %

    %
    % When running the code, determine a suitable time step. A too small time
    % step will make the calculation very long, while for a too large time step
    % the solution will blow up due to numerical instability.
    %

    Re = 1000;             % Reynolds number
%     N = 32;                % Number of volumes in the x- and y-direction
    Delta = 1/N;           % uniforme ruimte stap in de afbeelding

    % dt = 1.195^-N;           % time step

    a =    0.001747;
    b =  -2.441e-05;
    c =        1.13;
    d =      -0.1854;
    dt = a + b*N + c*exp(d*N);


%     tol = 1e-3;            % tol determines when steady state is reached and the program terminates

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


    u_diff = 1;
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

    while u_diff > tol
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

        u_diff = max(abs(u-uold))/dt; 
%         diff = 1e-7 % used in debugging 


%             figure(1)
%             plot(iter/1000, u_diff, 'bx')
%             hold on

        diff_list(end+1) = u_diff;


        if mod(iter,1000) == 0
            disp(u_diff)
            toc
        end
        

        iter = iter + 1;
        elapsed_time = toc;
    end
    

    % checks on pressure matrix from assignment

    if sum(A,2) < 1e-10
        disp('CHECK passed, row sum equals zero')
    end
    if A == transpose(A)
        disp('CHECK passed, A == transpose(A)')
    end

    %
    %  Produce desired out put to compare yur results with those given in the
    %  reference by Bptella and Peyret
    %

    % load Botella and Peyret values
    bot = load_bot();


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST-PRO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % make h_vector to convert circulation to velocity:
    h_vec = [h(i_u) h(j_v)]';
    vel_d1 = u./h_vec;                       % velocity along dual grid 1 forms 

    % dual grid 0-form i and j indices. 
    i_d0 = repmat(1:N, 1, N)';
    j_d0 = ceil((1:N*N)/N)';


    u_d0 = (findu(i_d0,j_d0, N, vel_d1, col_idx_tE21) + findu(i_d0+1, j_d0, N, vel_d1, col_idx_tE21))/2;
    v_d0 = (findv(i_d0,j_d0, N, vel_d1, col_idx_tE21) + findv(i_d0, j_d0+1, N, vel_d1, col_idx_tE21))/2;


    % p = P - 0.5*sqrt((u_d0.^2 + v_d0.^2)).^2; % ?? maybe one of these will give
    % the right pressure contours
    p = P - 0.5*(u_d0.^2 + v_d0.^2);

    % subtract reference pressure as in botella (mid point pressure)
    if mod(N,2)==0      % if N is even use the average of the four mid pressures
        p_mat = reshape(p,N,N)';
        mid_pressure = (p_mat(N/2,N/2)+p_mat(N/2,N/2+1)+p_mat(N/2+1,N/2)+p_mat(N/2+1,N/2+1))/4;
        p_mat = p_mat-mid_pressure;
    else                % if N is uneven use the mid pressure
        p = p - p(ceil(length(p)/2)); 
        p_mat = reshape(p,N,N)';
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VORTICITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    txi_mat = reshape(txi,N+1,N+1)'; % vorticity on primal grid 0-nodes in matrix 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STREAMLINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % convert circulation on dual grid to flux on primal grid: 
    tu = Ht11*u;

    % make matrix of u fluxes
    tu_mat = zeros(N,N+1);
    for ii = 1:N+1
        j_col = flip((1:N)');
        i_col = ii*ones(size(j_col));
        tu_mat(:,ii) = findu(i_col,j_col, N, tu, col_idx_tE21);
    end

    % % multiply fluxes by their respective dy (vertical edge lengths) already
    % done in u matrix?
    % tuh_mat = tu_mat.*repmat(th',1, N+1);

    % take cumulitive sum to represent integral
    psi_mat = cumsum(tu_mat, 'reverse');

    % add empty row to represent the bottom row nodes 
    psi_mat = flip([psi_mat;zeros(size(psi_mat(1,:)))]);


    %%%%%%%%%%%%% SAVE VERTICAL AND HORIZONTAL CENTRELINE RESULTS %%%%%%%%%%%%%

    u_d0_mat = reshape(u_d0,N,N)';
    v_d0_mat = reshape(v_d0,N,N)';


    sol.constx_y = x(2:end-1);
    sol.constx_y_xi = tx;

    sol.consty_x = x(2:end-1);
    sol.consty_x_xi = tx;

    if mod(N,2) == 0
        %vertical
        sol.constx_u = (u_d0_mat(:,N/2) + u_d0_mat(:,N/2+1))/2;
        sol.constx_p = (p_mat(:,N/2) + p_mat(:,N/2+1))/2;
        sol.constx_xi = txi_mat(:,N/2+1);
        %horizontal
        sol.consty_v = (v_d0_mat(N/2,:) + v_d0_mat(N/2+1, :))/2;
        sol.consty_p = (p_mat(N/2, :) + p_mat(N/2+1, :))/2;
        sol.consty_xi = txi_mat(N/2+1, :);    
    else
        %vertical
        sol.constx_u = u_d0_mat(:,ceil(N/2));
        sol.constx_p = p_mat(:,ceil(N/2));
        sol.constx_xi = (txi_mat(:,ceil(N/2)) + txi_mat(:,ceil(N/2)+1))/2;
        %horizontal
        sol.constx_v = v_d0_mat(ceil(N/2), :);
        sol.constx_p = p_mat(ceil(N/2),:);
        sol.constx_xi = (txi_mat(ceil(N/2),:) + txi_mat(ceil(N/2)+1,:))/2;    
    end


    totaldiff = 0;
    constx_u_interp = interp1(sol.constx_y, sol.constx_u, bot.constx_y());
    totaldiff = totaldiff + nansum(abs(bot.constx_u - constx_u_interp));

    constx_p_interp = interp1(sol.constx_y, sol.constx_p, bot.constx_y());
    totaldiff = totaldiff + nansum(abs(bot.constx_p - constx_p_interp));

    constx_xi_interp = interp1(sol.constx_y_xi, sol.constx_xi, bot.constx_y());
    totaldiff = totaldiff + nansum(abs(bot.constx_xi - constx_xi_interp));


    consty_v_interp = interp1(sol.consty_x, sol.consty_v, bot.consty_x());
    totaldiff = totaldiff + nansum(abs(bot.consty_v - consty_v_interp));

    consty_p_interp = interp1(sol.consty_x, sol.consty_p, bot.consty_x());
    totaldiff = totaldiff + nansum(abs(bot.consty_p - consty_p_interp));

    consty_xi_interp = interp1(sol.consty_x_xi, sol.consty_xi, bot.consty_x());
    totaldiff = totaldiff + nansum(abs(bot.consty_xi - consty_xi_interp));
end
