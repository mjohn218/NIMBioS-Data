function f = ktnumtobeintegrated(x,t,s,D,k)

h = 2.0 * pi * s * D;
alp = h * x .* j1(x * s) + k * j0(x * s);
bet = h * x .* y1(x * s) + k * y0(x * s);
tet = alp .* alp + bet .* bet;
P = j0(x * s) .* y1(x * s) - j1(x * s) .* y0(x* s);
T = (j0(x * s) .* bet - y0(x * s) .* alp) ./ tet;

f = T .* P .* (1 - exp(-D * t * x .* x));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function r = j0(inp)

r = besselj(0,inp);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function r = j1(inp)

r = besselj(1,inp);

function r = y0(inp)

r = bessely(0,inp);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function r = y1(inp)

r = bessely(1,inp);
