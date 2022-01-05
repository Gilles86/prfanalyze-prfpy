import six
import pimms
import multiprocessing
# import sharedmem
from matplotlib import pyplot as plt
from prfpy.timecourse import *
from prfpy.rf import *
from prfpy.fit import CFFitter, CSS_Iso2DGaussianFitter, DoG_Iso2DGaussianFitter, Iso2DGaussianFitter, Norm_Iso2DGaussianFitter
from prfpy.model import CFGaussianModel, CSS_Iso2DGaussianModel, DoG_Iso2DGaussianModel, Iso2DGaussianModel, Norm_Iso2DGaussianModel
from prfpy.stimulus import CFStimulus, PRFStimulus2D
# from bids.reports import BIDSReport
# from bids import BIDSLayout
import json
import pandas as pd
import nibabel as nib
import numpy as np
import os
import sys
from scipy.stats.stats import pearsonr

# Load and open files
opts_file = "/tank/mueller/projects/prfanalyze-prfpy/default_config.json"
bold_file = "/tank/mueller/projects/data/simulated/prfsynth_20_voxels/BIDS/sub-001/ses-20200320/func/sub-001_ses-20200320_task-prf_acq-normal_run-01_bold.nii.gz"
stim_file = "/tank/mueller/projects/data/simulated/prfsynth_20_voxels/BIDS/stimuli/sub-001_ses-20200320_task-prf_apertures.nii.gz"
stimjs_file = "/tank/mueller/projects/data/simulated/prfsynth_20_voxels/BIDS/derivatives/prfsynth/sub-001/ses-20200320/sub-001_ses-20200320_task-prf_acq-normal_run-01_bold.json"
outdir = "/tank/mueller/projects/data/simulated/prfsynth_20_voxels/BIDS/derivatives/prfanalyze-prfpy/sub-001/ses-20200320/"
# (opts_file, bold_file, stim_file, stimjs_file, outdir) = sys.argv[1:]

with open(opts_file, 'r') as fl:
    opts = json.load(fl)
    opts = opts.get('options')
    if opts is None:
        raise ValueError(
            'Please make sure the config json has the correct structure and contains an "options" section!')
bold_im = nib.load(bold_file)
stim_im = nib.load(stim_file)
with open(stimjs_file, 'r') as fl:
    stim_json = json.load(fl)


# seed random number generator so we get the same answers ...
# TODO how much randomness is used in PRFPY?
np.random.seed(opts.get('seed', 2764932))

# if there is an HRF search or not...
fixed_hrf = opts.get('fixed_hrf', False)


# Put data into correct shape
bold = np.reshape(np.asarray(bold_im.dataobj), (-1, bold_im.shape[-1]))
stim = np.squeeze(np.asarray(stim_im.dataobj))
if len(stim_json) != bold.shape[0]:
    raise ValueError(
        'BOLD Image and Stimulus JSON do not have the same number of data points')


if fixed_hrf:
    # fields = ('theta', 'rho', 'sigma', 'beta', 'baseline')
    fields = ('mu_x', 'mu_y', 'sigma', 'beta', 'baseline')
else:
    # fields = ('theta', 'rho', 'sigma', 'hrfdelay', 'beta', 'baseline')
    fields = ('mu_x', 'mu_y', 'sigma', 'beta', 'baseline', 'hrf_1', 'hrf_2')


class FittingConfig:
    """
    FittingConfig

    Class to configure everything needed to create the model and run the fitting
    """

    def __init__(self, opts: dict, fixed_hrf: bool, stim):
        # MULTIPROCESSING
        opt_multiprocessing = opts.get('multiprocess', False)
        if opt_multiprocessing == 'auto' or opt_multiprocessing:
            self.number_jobs = multiprocessing.cpu_count()
        elif type(opt_multiprocessing) is int and opt_multiprocessing > 1:
            self.number_jobs = opt_multiprocessing
        else:
            self.number_jobs = 1

        # STIMULUS
        self.dist = opts.get('screen_distance', 100)
        self.stim = stim

        # MODEL
        model_opts = opts.get('model', {})
        self.model_type = model_opts.get("model_type", "Iso2DGaussian")
        self.fit_hrf = not fixed_hrf

        # FITTING PARAMS
        # grid search
        fit_opts = opts.get('fitting', {})
        # TODO increase the default grid size
        self.grid_size = fit_opts.get('grid_size', 10)
        self.eccs_lower = fit_opts.get('eccs_lower',  0.01)
        self.eccs_upper = fit_opts.get('eccs_upper',  1.0)
        self.polars_lower = fit_opts.get('polars_lower',  5.0)
        self.polars_upper = fit_opts.get('polars_upper',  12.0)
        self.sizes_lower = fit_opts.get('sizes_lower',  0.01)
        self.sizes_upper = fit_opts.get('sizes_upper',  1.0)
        self.grid_size = fit_opts.get('grid_size',  10)
        # interative search
        self.eps = fit_opts.get('eps', 1e-1)
        self.xtol = fit_opts.get('xtol', 0.00001)
        self.ftol = fit_opts.get('ftol', 0.00001)
        self.constraints = []  # TODO check this
        self.rsq_threshold = fit_opts.get('rsq_threshold', 0.001)

        self._init_search_grid()

        # self.max_ecc_size = max(sizes)

        # TODO which bounds to use
        # gauss_bounds = [
        #     (-1.5 * max_ecc_size, 1.5 * max_ecc_size),  # mu_x
        #     (-1.5 * max_ecc_size, 1.5 * max_ecc_size),  # mu_y
        #     (eps, 1.5 * ss),                           # size
        #     (0, 1000),                                 # beta
        #     (0, 1000),                                 # baseline
        #     # (0, 1000),                               # hrf_1
        #     # (eps, 3*ss),                             # hrf_2
        # ]

    def init_stimulus_params(self, info: dict):
        self.stdat = info['Stimulus']
        if pimms.is_list(self.stdat):
            self.stdat = self.stdat[0]
        self.height = self.stdat['fieldofviewVert']
        self.width = self.stdat['fieldofviewHorz']
        self.stim_width = 2 * self.dist * np.tan(np.pi/180 * self.width/2)
        self.tr = info['TR']

    def _init_search_grid(self):
        # TODO Tomas: this is dependent on the stimulus but not on the different voxels right?
        self.eccs = np.linspace(
            self.eccs_lower, self.eccs_upper, self.grid_size, dtype='float32')
        self.polars = np.linspace(
            self.polars_lower, self.polars_upper, self.grid_size, dtype='float32')
        self.sizes = np.linspace(
            self.sizes_lower, self.sizes_upper, self.grid_size, dtype='float32')

    def get_search_grid(self):
        return self.eccs, self.polars, self.sizes

    def init_stimulus(self):
        if (self.model_type == "CFGaussian"):
            # CFStimulus(data=self.stim)
            pass
        # elif (self.model_type == "CSS_Iso2DGaussian"):

        # elif (self.model_type == "DoG_Iso2DGaussian"):
        # elif (self.model_type == "Norm_Iso2DGaussian"):
        else:  # (self.model_type == "Iso2DGaussian")
            self.stimulus = PRFStimulus2D(screen_size_cm=self.height,
                                          screen_distance_cm=self.dist,
                                          design_matrix=self.stim,
                                          TR=self.tr)
            #  task_lengths=task_lengths,
            #  task_names=task_names)

    def init_model(self):
        if (self.model_type == "CFGaussian"):
            # TODO implement CFGaussian
            self.gg = CFGaussianModel(stimulus=self.stimulus)
        elif (self.model_type == "CSS_Iso2DGaussian"):
            # TODO implement CSS_Iso2DGaussian
            self.gg = CSS_Iso2DGaussianModel(stimulus=self.stimulus)
        elif (self.model_type == "DoG_Iso2DGaussian"):
            # TODO implement DoG_Iso2DGaussian
            self.gg = DoG_Iso2DGaussianModel(stimulus=self.stimulus)
        elif (self.model_type == "Norm_Iso2DGaussian"):
            # TODO implement Norm_Iso2DGaussian
            self.gg = Norm_Iso2DGaussianModel(stimulus=self.stimulus)

        else:  # (self.model_type == "Iso2DGaussian")
            self.gg = Iso2DGaussianModel(stimulus=self.stimulus)

    def get_fitter(self, data: np.ndarray):
        if (self.model_type == "CFGaussian"):
            # TODO implement CFGaussian
            self.gf = CFFitter(data=data, model=self.gg,
                               fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        elif (self.model_type == "CSS_Iso2DGaussian"):
            # TODO implement CSS_Iso2DGaussian
            self.gf = CSS_Iso2DGaussianFitter(
                data=data, model=self.gg, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        elif (self.model_type == "DoG_Iso2DGaussian"):
            # TODO implement DoG_Iso2DGaussian
            self.gf = DoG_Iso2DGaussianFitter(
                data=data, model=self.gg, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        elif (self.model_type == "Norm_Iso2DGaussian"):
            # TODO implement Norm_Iso2DGaussian
            self.gf = Norm_Iso2DGaussianFitter(
                data=data, model=self.gg, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)

        else:  # (self.model_type == "Iso2DGaussian")
            self.gf = Iso2DGaussianFitter(
                data=data, model=self.gg, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
                
        return self.gf


def fit_voxel(x, info, config):
    config.init_stimulus_params(info)
    config.init_stimulus()
    config.init_model()

    # FITTER
    gf = config.get_fitter(data=x)

    # SEARCH GRID
    eccs, polars, sizes = config.get_search_grid()

    # Iterative fit
    max_ecc_size = config.stimulus.screen_size_degrees

    # TODO which bounds to use
    # TODO put the bound configuration also in the json file
    gauss_bounds = [
        (-1.5 * max_ecc_size, 1.5 * max_ecc_size),  # mu_x
        (-1.5 * max_ecc_size, 1.5 * max_ecc_size),  # mu_y
        # TODO is this correct?  # size
        (config.eps, 1.5 * config.stim_width),
        (0, 1000),                                 # beta
        (0, 1000),                                 # baseline
    ]

    if config.fit_hrf:
        gauss_bounds.append((
            (0, 1000),                               # hrf_1
            # hrf_2
            (config.eps, 3*config.stimulus.screen_size_degrees),
        ))

    # FIT
    gf.grid_fit(ecc_grid=eccs,
                polar_grid=polars,
                size_grid=sizes,
                verbose=False)

    print(f"GRID FIT finished successfully")

    gf.iterative_fit(rsq_threshold=config.rsq_threshold,
                     bounds=gauss_bounds,
                     constraints=config.constraints,
                     verbose=False)

    # Get RESULTS
    # params = gf.gridsearch_params[index, :]
    params = gf.iterative_search_params
    
    # print(f"Fit for vox {index} = {params[-1]}")
    pred = [config.gg.return_prediction(
        *params[i, :-1])[0] for i in range(len(params))]

    # this should return an aray with
    # index, voxel_data, estimates for all the fields, prediction
    return (list(range(len(x))), x) + tuple(params.T) + (pred,)


# Initialise the configuration containing the stimulus, the model and parameters for the fitting
fitting_config = FittingConfig(opts=opts, stim=stim, fixed_hrf=fixed_hrf)

voxs = fit_voxel(x=bold, info=stim_json[0], config=fitting_config)

# data_bundle = zip(range(len(bold)), bold, stim_json)
# if number_jobs == 1:
#     voxs = [fit_voxel((index, x, stim_json), fitting_config)
#             for (index, x, stim_json) in data_bundle]
# else:
#     print(f"Using {number_jobs} threads for parallel computing.")
#     tups = list(data_bundle)
#     with sharedmem.Pool(np=number_jobs) as pool:
#         voxs = pool.map(fit_voxel, [tups, fitting_config])
#     voxs = list(sorted(voxs, key=lambda tup: tup[0]))

# Update the results to match the x0/y0, sigma style used by prfanalyze
all_fields = ('index', 'voxel') + fields + ('fit', 'pred',)
res = {field_name: voxs[index]
       for (index, field_name) in enumerate(all_fields)}


r2s = []
for i in range(len(bold)):
    r2 = pearsonr(res['pred'][i], bold[i])
    print(r2)
    r2s.append(r2)

final_res = {}
# TODO since this is in degrees of visual angle, does it have to be converted back? If yes, to what?
final_res['centerx0'] = res['mu_x']
final_res['centery0'] = res['mu_y']
# final_res['centerx0'] = np.cos(res['theta']) * res['rho']
# final_res['centery0'] = -np.sin(res['theta']) * res['rho']
final_res['sigmamajor'] = res['sigma']
final_res['sigmaminor'] = res['sigma']
final_res['beta'] = res['beta']
final_res['baseline'] = res['baseline']


# Save results to files
for (attr_name, attr_data) in six.iteritems(final_res):
    im = nib.Nifti1Image(np.reshape(
        attr_data, bold_im.shape[:-1]), bold_im.affine)
    im.to_filename(os.path.join(outdir, attr_name + '.nii.gz'))
# Also export the prediction and testdata
im = nib.Nifti1Image(bold.reshape(bold_im.shape), bold_im.affine)
im.to_filename(os.path.join(outdir, 'testdata.nii.gz'))
im = nib.Nifti1Image(np.reshape(res['pred'], bold_im.shape), bold_im.affine)
im.to_filename(os.path.join(outdir, 'modelpred.nii.gz'))

# Store R2
# im = nib.Nifti1Image(r2s, bold_im.affine)
# im.to_filename(os.path.join(outdir, 'r2.nii.gz'))


# That's it!
print("PRFPY finished succesfully.")
sys.exit(0)
