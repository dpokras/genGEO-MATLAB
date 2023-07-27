function [C_surfacePipe] = CapitalCost_SurfacePipe(params)

PPI = PPI('PPI_Pipe', params.costYear); %2022
PPI_ref = PPI('PPI_Pipe', 2003)
X_PCs = 1.15;
X_ICs = 1.12;

if (strcmp(params.wellFieldType,'Doublet'))
    L_surfacePipe = 707;
    D_surfacePipe = 0.41;
elseif (strcmp(params.wellFieldType,'5spot'))
    L_surfacePipe = 3000;
    D_surfacePipe = 0.41;
elseif (strcmp(params.wellFieldType,'5spot_SharedNeighbor'))
    L_surfacePipe = 707;
    D_surfacePipe = 0.41;
elseif (strcmp(params.wellFieldType,'Tungsten'))
    L_surfacePipe = 800;
    D_surfacePipe = 0.41;
elseif (strcmp(params.wellFieldType,'5spot_ManyN'))
    switch params.N_5spotd
        case 1
            L_surfacePipe = 3000;
            D_surfacePipe = 0.41;
        case 2
            L_surfacePipe = 12000;
            D_surfacePipe = 0.54;
        case 3
            L_surfacePipe = 25000;
            D_surfacePipe = 0.65;
        case 4
            L_surfacePipe = 45000;
            D_surfacePipe = 0.79;
        case 5
            L_surfacePipe = 69000;
            D_surfacePipe = 0.89;
        case 6
            L_surfacePipe = 107000;
            D_surfacePipe = 0.98;
        case 7
            L_surfacePipe = 153000;
            D_surfacePipe = 1.02;
        case 8
            L_surfacePipe = 223000;
            D_surfacePipe = 1.05;
        case 9
            L_surfacePipe = 309000;
            D_surfacePipe = 1.09;
        case 10
            L_surfacePipe = 406000;
            D_surfacePipe = 1.13;
        otherwise
            throw(MException('CapitalCost_SurfacePipe:NotImplemented','Not Implemented'));
    end
else
    throw(MException('CapitalCost_SurfacePipe:NotImplemented','Not Implemented'));
end

c_surfacePipe = 2205 * D_surfacePipe^2 + 134;
C_surfacePipe = X_PCs * X_ICs * PPI/PPI_ref * c_surfacePipe * L_surfacePipe;

end

