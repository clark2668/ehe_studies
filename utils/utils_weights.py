import simweights
import nuflux

H3a_SIBYLL23C = ['H3a_SIBYLL23C', 'H3a_SIBYLL23C_pr', 'H3a_SIBYLL23C_conv',
				'H3a_SIBYLL23C_k', 'H3a_SIBYLL23C_K0', 'H3a_SIBYLL23C_K0L',
				'H3a_SIBYLL23C_K0S', 'H3a_SIBYLL23C_pi', 'H3a_SIBYLL23C_mu'
]

nuflux_models = H3a_SIBYLL23C

# a dictionary of custom fluxes
custom_nu_models = {}

def get_neutrino_flux_model(model):
	flux = None
	if model in nuflux_models:
		flux = nuflux.makeFlux(model)
	elif model in custom_nu_models.keys():
		flux = custom_nu_models[model]
	else:
		raise NotImplementedError('Flux model {} is not implemented'.format(model))
	return flux

cr_flux_models = {'GaisserH4a' : simweights.GaisserH4a(),
					'GaisserH3a' : simweights.GaisserH3a()
				}

def get_cr_flux_model(model):
	flux = cr_flux_models[model] # will raise KeyError, so should be transparent
	return flux

options = ['nugen', 'corsika']

def get_flux_model(model, simtype):
	if simtype is 'nugen':
		flux_model = get_neutrino_flux_model(model)
	elif simtype is 'corsika':
		flux_model = get_cr_flux_model(model)
	else:
		raise NotImplementedError('Simi type {} is not implemented'.format(simtype))
	return flux_model

def get_weighter(open_hdf_file, simtype, nfiles=1):
	'''
	please first open the file by doing

	open_hdf_file = pandas.HDFStore(file, 'r')

	'''
	if simtype is 'nugen':
		weighter = simweights.NuGenWeighter(open_hdf_file, nfiles=nfiles)
	elif simtype is 'corsika':
		weighter = simweights.CorsikaWeighter(open_hdf_file, nfiles=nfiles)
	else:
		raise NotImplementedError('Simi type {} is not implemented'.format(simtype))
	return weighter





