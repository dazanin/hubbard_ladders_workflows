from scripts import load, utils

all_sizes = [32, 48, 64, 80, 96, 128, 160, 192]
all_filling = [0.875, 0.9375, 0.96875]

all_corr_types = []
all_corr_types.append( utils.Averaged() )
all_corr_types.append( utils.FixedStart(18) )
all_corr_types.append( utils.FixedStart(20) )

all_extraps = []
all_extraps.append(800)
all_extraps.append(utils.Extrapolation('variance', deg=2, num_points=None))

if False:
    ## Cache pair correlations
    for filling in all_filling:
        for L in all_sizes:
            for corr_type in all_corr_types:
                for bond_dim in all_extraps:
                    load.result('pairfield', L=L, filling=filling, bond_dim=bond_dim, correlation_type=corr_type)

    ## Cache density-density correlations
    for filling in all_filling:
        for L in all_sizes:
            for corr_type in all_corr_types:
                for bond_dim in all_extraps:
                    load.result('densdens', L=L, filling=filling, bond_dim=bond_dim, correlation_type=corr_type)

    ## Cache density and density fit
    for filling in all_filling:
        for L in all_sizes:
            for bond_dim in all_extraps:
                load.result('density', L=L, filling=filling, bond_dim=bond_dim)

    ## Cache energy
    for filling in all_filling:
        for L in all_sizes:
            for bond_dim in all_extraps:
                if isinstance(bond_dim, utils.Extrapolation):
                    load.extrapolation('energy', L=L, filling=filling, bond_dim=bond_dim)


## Odd sizes
all_sizes = [33, 49, 65, 81, 97, 129]
all_filling = [0.875]

## Cache density and density fit
for filling in all_filling:
    for L in all_sizes:
        for bond_dim in all_extraps:
            load.result('density', L=L, filling=filling, bond_dim=bond_dim)

## Cache energy
for filling in all_filling:
    for L in all_sizes:
        for bond_dim in all_extraps:
            if isinstance(bond_dim, utils.Extrapolation):
                load.extrapolation('energy', L=L, filling=filling, bond_dim=bond_dim)
