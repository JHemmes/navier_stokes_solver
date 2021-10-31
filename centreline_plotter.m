clear all
close all

N16 = load('output/N16_centreline_var.mat');
N32 = load('output/N32_centreline_var.mat');
N48 = load('output/N48_centreline_var.mat');
N56 = load('output/N56_centreline_var.mat');
N64 = load('output/N64_centreline_var.mat');

bot = load_bot();

%%%%%%%%%%%%%%%%%%%%% PLOT ALONG VERTICAL CENTRELINE %%%%%%%%%%%%%%%%%%%%%%
%%
width = 12;
height = 10;

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
scatter(bot.constx_y, bot.constx_u, 'filled', 'linewidth', 1)
hold on 
plot(N16.sol.constx_y, N16.sol.constx_u, 'linewidth', 1)
plot(N32.sol.constx_y, N32.sol.constx_u, '--', 'linewidth', 1)
plot(N48.sol.constx_y, N48.sol.constx_u, ':', 'linewidth', 1)
plot(N56.sol.constx_y, N56.sol.constx_u, '-.', 'linewidth', 1)
plot(N64.sol.constx_y, N64.sol.constx_u, 'linewidth', 1)
legend('Benchmark results','$N = 16$','$N = 32$','$N = 48$','$N = 56$','$N = 64$','Location', 'SouthWest', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$y$', 'Interpreter', 'LaTeX')
ylabel('$u$', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/centre_vertical_u

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
scatter(bot.constx_y, bot.constx_p, 'filled', 'linewidth', 1)
hold on 
plot(N16.sol.constx_y, N16.sol.constx_p, 'linewidth', 1)
plot(N32.sol.constx_y, N32.sol.constx_p, '--', 'linewidth', 1)
plot(N48.sol.constx_y, N48.sol.constx_p, ':', 'linewidth', 1)
plot(N56.sol.constx_y, N56.sol.constx_p, '-.', 'linewidth', 1)
plot(N64.sol.constx_y, N64.sol.constx_p, 'linewidth', 1)
legend('Benchmark results','$N = 16$','$N = 32$','$N = 48$','$N = 56$','$N = 64$','Location', 'NorthEast', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$y$', 'Interpreter', 'LaTeX')
ylabel('$p$', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/centre_vertical_p

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
scatter(bot.constx_y, bot.constx_xi, 'filled', 'linewidth', 1)
hold on 
plot(N16.sol.constx_y_xi, N16.sol.constx_xi, 'linewidth', 1)
plot(N32.sol.constx_y_xi, N32.sol.constx_xi, '--', 'linewidth', 1)
plot(N48.sol.constx_y_xi, N48.sol.constx_xi, ':', 'linewidth', 1)
plot(N56.sol.constx_y_xi, N56.sol.constx_xi, '-.', 'linewidth', 1)
plot(N64.sol.constx_y_xi, N64.sol.constx_xi, 'linewidth', 1)
legend('Benchmark results','$N = 16$','$N = 32$','$N = 48$','$N = 56$','$N = 64$','Location', 'NorthWest', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$y$', 'Interpreter', 'LaTeX')
ylabel('$\xi$', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/centre_vertical_xi


% %%%%%%%%%%%%%%%%%%%%% PLOT ALONG HORIZONTAL CENTRELINE %%%%%%%%%%%%%%%%%%%%

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
scatter(bot.consty_x, bot.consty_v, 'filled', 'linewidth', 1)
hold on 
plot(N16.sol.consty_x, N16.sol.consty_v, 'linewidth', 1)
plot(N32.sol.consty_x, N32.sol.consty_v, '--', 'linewidth', 1)
plot(N48.sol.consty_x, N48.sol.consty_v, ':', 'linewidth', 1)
plot(N56.sol.consty_x, N56.sol.consty_v, '-.', 'linewidth', 1)
plot(N64.sol.consty_x, N64.sol.consty_v, 'linewidth', 1)
legend('Benchmark results','$N = 16$','$N = 32$','$N = 48$','$N = 56$','$N = 64$','Location', 'SouthEast', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$x$', 'Interpreter', 'LaTeX')
ylabel('$v$', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/centre_horizontal_v

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
scatter(bot.consty_x, bot.consty_p, 'filled', 'linewidth', 1)
hold on 
plot(N16.sol.consty_x, N16.sol.consty_p, 'linewidth', 1)
plot(N32.sol.consty_x, N32.sol.consty_p, '--', 'linewidth', 1)
plot(N48.sol.consty_x, N48.sol.consty_p, ':', 'linewidth', 1)
plot(N56.sol.consty_x, N56.sol.consty_p, '-.', 'linewidth', 1)
plot(N64.sol.consty_x, N64.sol.consty_p, 'linewidth', 1)
legend('Benchmark results','$N = 16$','$N = 32$','$N = 48$','$N = 56$','$N = 64$','Location', 'North', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$x$', 'Interpreter', 'LaTeX')
ylabel('$p$', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/centre_horizontal_p

figure('Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
scatter(bot.consty_x, bot.consty_xi, 'filled', 'linewidth', 1)
hold on 
plot(N16.sol.consty_x_xi, N16.sol.consty_xi, 'linewidth', 1)
plot(N32.sol.consty_x_xi, N32.sol.consty_xi, '--', 'linewidth', 1)
plot(N48.sol.consty_x_xi, N48.sol.consty_xi, ':', 'linewidth', 1)
plot(N56.sol.consty_x_xi, N56.sol.consty_xi, '-.', 'linewidth', 1)
plot(N64.sol.consty_x_xi, N64.sol.consty_xi, 'linewidth', 1)
legend('Benchmark results','$N = 16$','$N = 32$','$N = 48$','$N = 56$','$N = 64$','Location', 'South', 'Interpreter', 'LaTeX')
grid on 
grid minor
xlabel('$x$', 'Interpreter', 'LaTeX')
ylabel('$\xi$', 'Interpreter', 'LaTeX');
set(gca,'FontSize',12)
print -depsc2 figures/centre_horizontal_xi












