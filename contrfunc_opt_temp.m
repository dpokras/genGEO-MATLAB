function [c,ceq] = contrfunc(result, params)
c(1) = 1.6e3 - result.s_turb_in;
c(2) = params.side_stream_radius - params.well_radius;
ceq = result.W_net / 1e3 - 1.2e3;
end