function [c,ceq] = contrfunc(result)
c = 1.65e3 - result.s_turb_in;
% c = 100 - result.T_prod_surface;
ceq = [];
% ceq = result.W_net / 1e3 - 1.2e3;
end