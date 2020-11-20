function r = hyperbola_critical(x)
% compute the (real) critical points of dist(x, H), where H = {(s,t) : st=1}.
% returns a vector of at most 4 points r(i) s.t. (r(i), 1/r(i)) are these
% critical points.

r = roots([1 -x(1) 0 x(2) -1]);

r = r(imag(r) == 0); %excludes complex solutions

assert(length(r) >= 2); % there should be at least 2 real roots to this equation, by Descartes' rule of signs.
