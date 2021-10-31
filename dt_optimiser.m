

%% This section finds the convergence times for different timesteps 
dt_list = 0.001:0.001:0.1;
N = 16;

timelist = zeros(1,length(dt_list));
for ii = 1:length(dt_list)
    disp(dt_list(ii))
    timelist(ii) = funNS_solver_dt(dt_list(ii), N);
end

figure()
plot(dt_list(timelist ~= 10000), timelist(timelist ~= 10000), 'linewidth', 1)
grid on
% legend('$A(\eta)$','$B(\eta)$','Location', 'SouthWest', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$\Delta t$ [s]', 'Interpreter', 'LaTeX')
ylabel('Convergence time [s]', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/ConvergenceTime




%% This section is used to find the maximum dt for convergence

dt_list = flip(5e-4:1e-5:5e-3); % change this to values close to the expected dt
N = 32;

% note the loop can be manually aborted when the program starts converging
% to avoid having to wait for the solution. 
for ii = 1:length(dt_list)
    disp(dt_list(ii))
    duration = funNS_solver_dt(dt_list(ii), N);
    disp(duration)
    if duration < 10000
        duration
        dt_list(ii)
        break 
    end
end


%% This section fits the curve to the N-dt value pairs found above. 

dtlist=[0.05959 0.0040 0.00079 0.00042 0.00025]';
nlist=[16 32 48 56 64]';

g = fittype('a + b*x + c*exp(d*x)');
f0 = fit(nlist,dtlist,g, 'start', [0 1 3 -0.2]);

a = f0.a - 4e-5; % 4e-5 correction makes sure dt lies below the maximum allowable dt 
b = f0.b;
c = f0.c;
d = f0.d;

nlarge = 16:1:64;
dtlarge = a + b*nlarge + c*exp(d*nlarge);


figure()
scatter(nlist, dtlist, 'filled')
hold on 
plot(nlarge, dtlarge, 'linewidth', 1)
grid on
legend('max $\Delta t$ for convergence','Fitted curve','Location', 'NorthEast', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$N$', 'Interpreter', 'LaTeX')
ylabel('$\Delta t$ [s]', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/maxdt

