if(model_SPE)
%SPE model
[nx, ny, nz] = deal(60, 220, numel(layers));
cartDims = [nx, ny,nz];
physDims = cartDims .* [20, 20, 2]*ft;
G = computeGeometry(cartGrid(cartDims, physDims));
% Construct the model
rock     = SPE10_rock(layers);
else
%% Layered model


%Create the grid
[nx, ny, nz] = deal(125, 125, 1);
cartDims = [nx, ny,nz];
physDims = cartDims .* [20, 20, 2]*ft;
%physDims = [nx, ny,nz];
G = computeGeometry(cartGrid(cartDims, physDims));

% Set up uniform permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*1;
rock.poro = ones(G.cells.num, 1)*0.2;

% Create layers of permeability
nlay = nx/5;
v = [];
for i = nlay+1:2*nlay:ny
    for j = 0:nlay-1
        if i+j < ny+1
            v = [v j+i];
        end
    end
end


[I] = Sub2ind_g([1:nx],v,1:nz,nx,ny,nz);
rock.perm(I) = 1*10^(per);
end

N = nx*ny*nz;
rock.perm = convertFrom(rock.perm, milli*darcy);
is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));