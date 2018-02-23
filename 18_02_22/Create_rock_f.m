function [cartDims, physDims, G, rock] = Create_rock_f(model_SPE,nx, ny, nz, Lx, Ly, Lz, varargin)
opt = struct(  'layers'    , 1   , ...
    'per'       , 0   , ...
    'szl'       , 5   , ...
    'l_dir'     , 'x');
opt      = merge_options(opt, varargin{:});
layers   = opt.layers;
per      = opt.per;
szl      = opt.szl;
l_dir    = opt.l_dir;
if(model_SPE)
    %SPE model
    cartDims = [nx, ny,nz];
    physDims = cartDims .* [Lx, Ly, Lz]*ft;
    G = computeGeometry(cartGrid(cartDims, physDims));
    % Construct the model
    rock     = SPE10_rock(layers);
    is_pos             = rock.poro > 0;
   rock.poro(~is_pos) = min(rock.poro(is_pos));
else
    %% Layered model
    %Create the grid
    %[nx, ny, nz] = deal(125, 125, 1);
    cartDims = [nx, ny,nz];
    physDims = cartDims .* [Lx, Ly, Lz]*ft;
    %physDims = [nx, ny,nz];
    G = computeGeometry(cartGrid(cartDims, physDims));
    
    % Set up uniform permeability and constant porosity
    rock.perm = ones(G.cells.num, 1)*10;
    rock.poro = ones(G.cells.num, 1)*0.2;
    
    % Create layers of permeability
    nlay = nx/szl;
    v = [];
    for i = nlay+1:2*nlay:ny
        for j = 0:nlay-1
            if i+j < ny+1
                v = [v j+i];
            end
        end
    end
    if(l_dir=='x')
        [I] = round(Sub2ind_g(1:nx,v,1:nz,nx,ny,nz));
    else if (l_dir=='y')
            [I] = round(Sub2ind_g(v,1:ny,1:nz,nx,ny,nz));
        else
            [I] = round(Sub2ind_g(1:nx,1:ny,v,nx,ny,nz));
        end
    end

    rock.perm(I) = 10*10^(-per);
    
    rock.perm = convertFrom(rock.perm, milli*darcy);
    is_pos             = rock.poro > 0;
    rock.poro(~is_pos) = min(rock.poro(is_pos));
end