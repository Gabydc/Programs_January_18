 G   = computeGeometry(cartGrid([3, 3, 5]));

   f   = initSingleFluid('mu' , 1*centi*poise, ...
                         'rho', 1000*kilogram/meter^3);
                     disp(f)
   rock.perm = rand(G.cells.num, 1)*darcy()/100;

   bc  = pside([], G, 'LEFT', 2*barsa);
   src = addSource([], 1, 1);
   W   = verticalWell([], G, rock, 1, G.cartDims(2), []   , ...
                      'Type', 'rate', 'Val', 1*meter^3/day, ...
                      'InnerProduct', 'ip_tpf');
   W   = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), [], ...
                      'Type', 'bhp', 'Val', 1*barsa, ...
                      'InnerProduct', 'ip_tpf');

   T   = computeTrans(G, rock);

   state         = initResSol (G, 10);
   state.wellSol = initWellSol(G, 10);

   state = incompTPFA(state, G, T, f, 'bc', bc, 'src', src, ...
                      'wells', W, 'MatrixOutput', true);

   plotCellData(G, state.pressure)