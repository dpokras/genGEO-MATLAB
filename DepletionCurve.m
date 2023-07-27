function [gamma] = DepletionCurve(psi, p1, p2, p3)

if (psi <= (-p1/p2))
    gamma = 1 - (1 - p3) * (1 - 0.5*erfc(p2*psi + p1));
else
    C = 1.13;
    gamma = 1 - (1 - p3) * (1 - 0.5*exp(-C*(p2*psi + p1)));
end

end

