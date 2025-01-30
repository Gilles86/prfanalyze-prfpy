import prfpy
import six
import pimms
import multiprocessing
from matplotlib import pyplot as plt
from prfpy.timecourse import *
from prfpy.rf import *
from prfpy.fit import CFFitter, CSS_Iso2DGaussianFitter, DoG_Iso2DGaussianFitter, Fitter, Iso2DGaussianFitter, Norm_Iso2DGaussianFitter
from prfpy.model import CFGaussianModel, CSS_Iso2DGaussianModel, DoG_Iso2DGaussianModel, Iso2DGaussianModel, Norm_Iso2DGaussianModel
from prfpy.stimulus import CFStimulus, PRFStimulus2D
if hasattr(prfpy.utils, 'HRF'):
    from prfpy.utils import HRF
import json
import ruamel.yaml as yaml
import pandas as pd
import nibabel as nib
import numpy as np
import os
import sys
from scipy.stats.stats import pearsonr


class FittingConfig:
    """Class to configure everything needed to create the model and run the fitting



        Example
        ----------
        >>> config.init_stimulus_params(info)
        >>> config.init_stimulus()
        >>> config.init_gauss_bounds()
        >>> config.init_model()
        >>> fitter = config.get_fitter(data=x)
        >>> grid_fit_params, iterative_fit_params = config.get_fitting_params()
        >>> fitter.grid_fit(**grid_fit_params)
        >>> fitter.iterative_fit(**iterative_fit_params)
        >>> params = fitter.iterative_search_params
        >>> pred = [config.model.return_prediction(
            *params[i, :-1])[0] for i in range(len(params))]
    """

    def __init__(self, opts: dict, stim: np.ndarray, output_dir: str, stimjs_file: str = ''):
        """__init__

            Initialize the configuration object for the fitting.

            Parameters
            ----------
            opts: dict
                Dictionary containing the options section from the configuration yaml file.
            stim: ndarray
                Array containing the stimulus data. Should be loaded from a nifti file, and should have 3 dimensions with shape (n, m, t) 
                where n and m correspond to horizontal and vertical screen size, respectively and t corresponds to the number of timesteps.

            Example
            ----------
            >>> fitting_config = FittingConfig(opts=opts, stim=stimulus) 
        """
        # OUTPUT DIRECTORY
        self.output_dir = output_dir

        # MULTIPROCESSING
        opt_multiprocessing = opts.get('multiprocess', False)
        if opt_multiprocessing == 'auto' or opt_multiprocessing == True:
            self.number_jobs = multiprocessing.cpu_count()
        elif type(opt_multiprocessing) is int and opt_multiprocessing > 1:
            self.number_jobs = opt_multiprocessing
        else:
            self.number_jobs = 1

        # VERBOSE
        self.verbose = opts.get('verbose', False)

        # STIMULUS
        self.dist = opts.get('screen_distance', 57) # at 57 cm distance, 1cm equals 1 degree of visual angle
        self.stim = stim

        # TODO change default to prfpy
        synth = opts['synth']
        if synth['origin'] == 'prfsynth' and len(stimjs_file) > 0:
            with open(stimjs_file, 'r') as fl:
                stim_json = json.load(fl)
            self.init_stimulus_params(stim_json)

        else:  # synth == 'prfpy'
            self.tr = synth.get('tr', 1)
            self.width = synth.get('width', 20)
            self.height = synth.get('height', 20)
            self.stimulus_width_degrees = synth.get(
                'stimulus_width_degrees', 2)

        # MODEL
        model_opts = opts.get('model', {})
        self.hrf_opts = model_opts.get('hrf', None)
        self.model_type = model_opts.get("model_type", "Iso2DGaussian")
        fixed_hrf = opts.get('fixed_hrf', True)
        self.fit_hrf = not fixed_hrf

        # TODO this will be dependent on the model
        if fixed_hrf:
            self.fields = ('mu_x', 'mu_y', 'sigma', 'beta', 'baseline')
        else:
            self.fields = ('mu_x', 'mu_y', 'sigma', 'beta',
                           'baseline', 'hrf_1', 'hrf_2')

        # remove this block when migrating to the new HRF model
        # if not hasattr(prfpy.utils, 'HRF'):
        #     if fixed_hrf:
        #         # TODO implement computation of hrf function given parameters
        #         if self.hrf_opts == 'vista_twogammas':
        #             self.hrf = self._get_default_vista_twogammas_hrf()
        #         else:
        #             self.hrf = None
        #     else:
        #         self.hrf = None

        # FITTING PARAMS
        # grid search
        fit_opts = opts.get('fitting', {})
        self.fit_opts = fit_opts
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

        # self._init_search_grid()

    def _get_default_vista_twogammas_hrf(self, params=[5.4, 5.2, 10.8, 7.35, 0.35], t=list(range(0, 20))):
        """ _get_default_vista_twogammas_hrf

            Return the computed HRF timecourse given by the vistasoft implementation. 
            DO NOT USE, this is not a verified method and does not ensure the validity of the computed HRF.

            Parameters
            ----------
            params: list, optional
                List of parameter in the following order: [peak1, fwhm1, peak2, fwhm2, dip], by default [5.4, 5.2, 10.8, 7.35, 0.35]
            t: list, optional
                List of timesteps., by default list(range(0,20))

            Returns
            ----------
            list
                Computed HRF using two gamma functions for the given parameters and timesteps

            Example
            ----------
            >>> self._get_default_vista_twogammas_hrf() 
        """
        #  params
        peak1 = params[0]
        fwhm1 = params[1]
        peak2 = params[2]
        fwhm2 = params[3]
        dip = params[4]

        #  sanity check
        if(peak1 == 0 or fwhm1 == 0):
            # self.hrf = None
            return None
        else:
            #  Taylor:
            alpha1 = np.power(peak1, 2)/np.power(fwhm1, 2)*8*np.log(2)
            beta1 = np.power(fwhm1, 2)/peak1/8/np.log(2)
            gamma1 = np.multiply(np.power(np.divide(t, peak1), alpha1), np.exp(
                np.divide(-1 * (np.subtract(t, peak1)), beta1)))

            if peak2 > 0 and fwhm2 > 0:
                alpha2 = np.power(peak2, 2)/np.power(fwhm2, 2)*8*np.log(2)
                beta2 = np.power(fwhm2, 2)/peak2/8/np.log(2)
                gamma2 = np.multiply(np.power(np.divide(t, peak2), alpha2), np.exp(
                    np.divide(-1 * (np.subtract(t, peak2)), beta2)))
            else:
                gamma2 = min(abs(np.subtract(t, peak2))) == abs(
                    np.subtract(t, peak2))

            h = np.subtract(gamma1, dip*gamma2)
            hrf = h.reshape((1, h.shape[0]))
            return hrf

    def init_gauss_bounds(self):
        """init_gauss_bounds

            Initialize the gauss_bounds needed for the iterative fit for the Iso2DGaussian model (or possibly others).
            Stimulus needs to be initialize first by calling init_stimulus().

            Example
            ----------
            >>> fitting_config.init_gauss_bounds() 
        """
        # check whether stimulus has been initilized
        assert self.stimulus != None

        # Iterative fit
        gauss_bounds = self.fit_opts.get('gauss_bounds', {})

        if self.verbose:
            print(f"Screen_distance_cm: {self.stimulus.screen_distance_cm}")
            print(
                f"Screen_size_degrees/max_ecc_size: {self.stimulus.screen_size_degrees}")
            print(f"Stim width: {self.stim_width}")

        max_ecc_size = self.stimulus.max_ecc
        ss = self.stimulus.screen_size_degrees

        # TODO which bounds to use
        self.gauss_bounds = [
            (gauss_bounds['mu_x']['lower_factor'] * ss,
             gauss_bounds['mu_x']['upper_factor'] * ss),  # mu_x
            (gauss_bounds['mu_y']['lower_factor'] * ss,
             gauss_bounds['mu_y']['upper_factor'] * ss),  # mu_y
            (gauss_bounds['size']['lower_factor'] * self.eps,
             #  gauss_bounds['size']['upper_factor'] * self.stim_width),
             gauss_bounds['size']['upper_factor'] * ss),
            (gauss_bounds['beta']['lower'], gauss_bounds['beta']
             ['upper']),                                 # beta
            (gauss_bounds['baseline']['lower'], gauss_bounds['baseline']
             ['upper']),                                 # baseline
        ]

        if self.fit_hrf:
            self.gauss_bounds.append(
                (gauss_bounds['hrf_1']['lower'], gauss_bounds['hrf_1']['upper']))                               # hrf_1
            # hrf_2
            self.gauss_bounds.append((gauss_bounds['hrf_2']['lower_factor'] * self.eps, gauss_bounds['hrf_2']['upper_factor'] * ss),
                                     )

        if self.verbose:
            print(f"Using gauss bounds:\n{self.gauss_bounds}\n")

    def init_hrf(self):
        # assert self.tr

        # if hasattr(prfpy.utils, 'HRF'):
        #     self.hrf = HRF()

        if not self.fit_hrf:
            # TODO implement computation of hrf function given parameters
            if self.hrf_opts == 'vista_twogammas':
                self.hrf.values = self._get_default_vista_twogammas_hrf()
            else:
                self.hrf.create_spm_hrf(TR=self.tr, force=True)
        else:
            self.hrf.create_spm_hrf(TR=self.tr, force=True)

    def init_stimulus_params(self, info):
        """init_stimulus_params

            Initialize the parameters needed to create the stimulus object. 

            Parameters
            ----------
            info: dict
                Dictionary providing the data for the stimulus parameters.

            Example
            ----------
            >>> fitting_config.init_stimulus_params(info) 
        """
        if pimms.is_list(info):
            info = info[0]
        self.stdat = info['Stimulus']
        if pimms.is_list(self.stdat):
            self.stdat = self.stdat[0]
        self.height = self.stdat['fieldofviewVert']
        self.width = self.stdat['fieldofviewHorz']
        # self.stim_width = 2 * self.dist * np.tan(np.pi/180 * self.width/2)
        self.stim_width = self.stdat['barWidth']
        self.tr = info['TR']

    def init_search_grid(self):
        """_init_search_grid

            Initialize the search grid given the lower and upper bounds for eccentricity, 
            polar angle and size and a number of points in the grid grid_size.
            Will directly be called at the end of the construction.

            Example
            ----------
            >>> self._init_search_grid()
        """
        # TODO Tomas: this is dependent on the stimulus but not on the different voxels right?
        max_ecc_size = self.stimulus.max_ecc
        self.eccs = max_ecc_size * np.linspace(
            self.eccs_lower, self.eccs_upper, self.grid_size, dtype='float32') ** 2
        self.polars = np.linspace(
            self.polars_lower, np.pi*self.polars_upper, self.grid_size, dtype='float32')
        self.sizes = max_ecc_size * np.linspace(
            self.sizes_lower, self.sizes_upper, self.grid_size, dtype='float32') ** 2
        # self.eccs = np.linspace(
        #     self.eccs_lower, self.eccs_upper, self.grid_size, dtype='float32')
        # self.polars = np.linspace(
        #     self.polars_lower, self.polars_upper, self.grid_size, dtype='float32')
        # self.sizes = np.linspace(
        #     self.sizes_lower, self.sizes_upper, self.grid_size, dtype='float32')

    def get_search_grid(self):
        """get_search_grid

            Return the computed search grid. If this has not been initialize, call _init_search_grid first.

            Returns
            ----------
            tuple
                Search grid for eccentricity, polar angle and size.

            Example
            ----------
            >>> eccs, polars, sizes = fitting_config.get_search_grid() 
        """
        if self.eccs is None or self.polars is None or self.sizes is None:
            self._init_search_grid()

        return self.eccs, self.polars, self.sizes

    def init_stimulus(self):
        """init_stimulus

            Initialize the stimulus objects depending on the used model. Default is the PRFStimulus2D object.
            Before calling this init_stimulus_params() has to be called to ensure that necessary information is present.

            Example
            ----------
            >>> fitting_config.init_stimulus()
        """
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

    def init_model(self):
        """init_model

            Initialize the model object given the specific model type.
            Possible models are CFGaussianModel, CSS_Iso2DGaussianModel, DoG_Iso2DGaussianModel, Norm_Iso2DGaussianModel, Iso2DGaussianModel.
            Model will be store in self.model

            Example
            ----------
            >>> fitting_config.init_model()
        """

        # check whether stimulus has been initialized
        assert self.stimulus != None

        if (self.model_type == "CFGaussian"):
            self.model = CFGaussianModel(stimulus=self.stimulus)
        elif (self.model_type == "CSS_Iso2DGaussian"):
            self.model = CSS_Iso2DGaussianModel(
                stimulus=self.stimulus, hrf=self.hrf)
        elif (self.model_type == "DoG_Iso2DGaussian"):
            self.model = DoG_Iso2DGaussianModel(
                stimulus=self.stimulus, hrf=self.hrf)
        elif (self.model_type == "Norm_Iso2DGaussian"):
            self.model = Norm_Iso2DGaussianModel(
                stimulus=self.stimulus, hrf=self.hrf)
        else:  # (self.model_type == "Iso2DGaussian")
            self.model = Iso2DGaussianModel(
                stimulus=self.stimulus, hrf=self.hrf)

    def get_fitter(self, data: np.ndarray) -> "Fitter":
        """get_fitter

            Initialize and return the fitter with the given data depending on the model type.
            Will perform parallel fitting depending on self.number_jobs.
            Will perform HRF fitting depending on self.fit_hrf.
            Possible fitter classes are:
            CFFitter
            CSS_Iso2DGaussianFitter
            DoG_Iso2DGaussianFitter
            Norm_Iso2DGaussianFitter
            Iso2DGaussianFitter

            Parameters
            ----------
            data: np.ndarray
                Array containing the data that should be used for the fitting.

            Returns
            ----------
            Fitter
                Fitting object specific to the model type and data.

            Example
            ----------
            >>> 
        """
        if (self.model_type == "CFGaussian"):
            # TODO implement CFGaussian
            self.fitter = CFFitter(data=data, model=self.model,
                                   fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        elif (self.model_type == "CSS_Iso2DGaussian"):
            # TODO implement CSS_Iso2DGaussian
            self.fitter = CSS_Iso2DGaussianFitter(
                data=data, model=self.model, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        elif (self.model_type == "DoG_Iso2DGaussian"):
            # TODO implement DoG_Iso2DGaussian
            self.fitter = DoG_Iso2DGaussianFitter(
                data=data, model=self.model, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        elif (self.model_type == "Norm_Iso2DGaussian"):
            # TODO implement Norm_Iso2DGaussian
            self.fitter = Norm_Iso2DGaussianFitter(
                data=data, model=self.model, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        else:  # (self.model_type == "Iso2DGaussian")
            self.fitter = Iso2DGaussianFitter(
                data=data, model=self.model, fit_hrf=self.fit_hrf, n_jobs=self.number_jobs)
        return self.fitter

    def get_fitting_params(self) -> "tuple[dict, dict]":
        """Evaluate the configuration and return two dictionaries containing the needed parameters to call grid_fit and iterative_fit.

            Returns
            ----------
            dict, dict
                Dictionaries containing parameters for grid_fit and iterative_fit respectively for the configured model type.

            Example
            ----------
            >>> grid_fit_params, iterative_fit_params = get_fitting_params()
        """
        if self.model_type == 'Iso2DGaussian':
            eccs, polars, sizes = self.get_search_grid()

            grid_fit_params = {"ecc_grid": eccs,
                               "polar_grid": polars,
                               "size_grid": sizes,
                               "verbose": self.verbose}

            iterative_fit_params = {"rsq_threshold": self.rsq_threshold,
                                    "bounds": self.gauss_bounds,
                                    "xtol": self.xtol,
                                    "ftol": self.ftol,
                                    "constraints": self.constraints,
                                    "verbose": self.verbose}
        # elif self.model_type == 'Norm_Iso2DGaussian':
        #     grid_fit_params = {"surround_amplitude_grid": "",
        #                        "surround_size_grid": "",
        #                        "neural_baseline_grid": "",
        #                        "surround_baseline_grid": "",
        #                        "verbose": self.verbose}
        else:
            grid_fit_params = None
            iterative_fit_params = None
            print(f"Only Iso2DGaussion model and fitter are so far implemented.")

        return grid_fit_params, iterative_fit_params

    def save_as_nifti(self, data, shape, affine, filename: str):
        """Save array like data to nifti image

            

            Parameters
            ----------
            data: array_like
                Array like data to be saved to nifti file.
            shape: shape
                Shape to reshape the data to.
            affine: None or (4,4) array_like
                Affine for Nifti1Image
            filename: str
                Name of the nifti file including file ending.

            Example
            ----------
            >>> fitting_config.save_to_nifit(data=data, shape=shape, affine=affine, filename=filename)
        """
        im = nib.Nifti1Image(np.reshape(data, shape), affine)
        im.to_filename(os.path.join(self.output_dir, filename))

    def save_final_results(self, results: dict, shape, affine):
        """Save the relevant recovered parameters (see self.fields) to individual nifti files.

            

            Parameters
            ----------
            results: dict
                Dictionary containing the recovered parameters as returned by the fitter.
            shape: shape
                Shape to reshape the per voxel parameter array_like data structures to.
            affine: None or (4,4) array_like
                Affine for the nifti1Image

            Example
            ----------
            >>> fitting_config.save_final_results(results=results, shape=shape, affine=affine)
        """
        final_res = {}

        if (self.model_type == "CFGaussian"):
            pass
        elif (self.model_type == "CSS_Iso2DGaussian"):
            pass
        elif (self.model_type == "DoG_Iso2DGaussian"):
            pass
        elif (self.model_type == "Norm_Iso2DGaussian"):
            pass
        else:  # (self.model_type == "Iso2DGaussian")
            final_res['centerx0'] = results['mu_x']
            final_res['centery0'] = results['mu_y']
            final_res['sigmamajor'] = results['sigma']
            final_res['sigmaminor'] = results['sigma']
            final_res['beta'] = results['beta']
            final_res['baseline'] = results['baseline']

        # Save results to files
        for (attr_name, attr_data) in six.iteritems(final_res):
            self.save_as_nifti(attr_data, shape, affine, f'{attr_name}.nii.gz')
            # im = nib.Nifti1Image(np.reshape(
            #     attr_data, shape), affine)
            # im.to_filename(os.path.join(outdir, attr_name + '.nii.gz'))

    def init_fitting(self):
        """Call all necessary submethods from this class 
        giving the initial configuration to set up the hrf, 
        stimulus, search grid, gauss bounds and model. 
        Can also be manually for a more customised fitting.

            

            Example
            ----------
            >>> fitting_config.init_fitting()
        """
        self.init_hrf()
        self.init_stimulus()
        self.init_search_grid()
        self.init_gauss_bounds()
        self.init_model()


def fit_voxels(x, config: FittingConfig):
    """Fit all present voxels using the given fitting configuration and information for each voxel.

        Parameters
        ----------
        x: ndarray
            Array containing the data for each voxel.
        config: FittingConfig
            FittingConfig object containing all parameters needed for the fitting

        Returns
        ----------
        List
            List containing indices, the original data, the fitted parameters and the resulting predictions using those parameters.

        Example
        ----------
        >>> voxels = fit_voxels(x=bold, info=info, config=fitting_config)
    """
    # Initialise everything left to start the prfpy fitting
    config.init_fitting()

    # FITTER
    fitter = config.get_fitter(data=x)

    # SEARCH GRID
    grid_fit_params, iterative_fit_params = config.get_fitting_params()

    assert grid_fit_params is not None and iterative_fit_params is not None

    if config.verbose:
        print(f"Starting GRID FIT")

    # FIT
    # pass items from dictionary as if they were method arguments
    fitter.grid_fit(**grid_fit_params)

    if config.verbose:
        print(f"GRID FIT finished successfully")
        print(f"Starting ITERATIVE FIT")

    fitter.iterative_fit(**iterative_fit_params)

    if config.verbose:
        print(f"ITERATIVE FIT finished successfully")

    # Get RESULTS
    params = fitter.iterative_search_params

    pred = [config.model.return_prediction(
        *params[i, :-1])[0] for i in range(len(params))]

    # this should return an aray with
    # index, voxel_data, estimates for all the fields, prediction
    results = {
        'indices': list(range(len(x))),
        'voxels': x
    }
    for index, field in enumerate(config.fields):
        results[field] = (params.T)[index]

    results['pred'] = pred
    # return (list(range(len(x))), x) + tuple(params.T) + (pred,)
    return results


if __name__ == "__main__":
    # Load and open files
    (opts_file, bold_file, stim_file, stimjs_file, outdir) = sys.argv[1:]

    with open(opts_file, 'r') as stream:
        try:
            yamls = yaml.YAML(typ='safe', pure=True)
            opts = yamls.load(stream)
        except yaml.YAMLError:
            raise ValueError(
                'Please make sure the config yaml has the correct structure!')

    bold_im = nib.load(bold_file)
    stim_im = nib.load(stim_file)

    # Put data into correct shape
    bold = np.reshape(np.asarray(bold_im.dataobj), (-1, bold_im.shape[-1]))
    stim = np.squeeze(np.asarray(stim_im.dataobj))

    # seed random number generator so we get the same answers ...
    # TODO how much randomness is used in PRFPY?
    np.random.seed(opts.get('seed', 2764932))

    # Initialise the configuration containing the stimulus, the model and parameters for the fitting
    fitting_config = FittingConfig(opts=opts, stim=stim, output_dir=outdir, stimjs_file=stimjs_file)

    # Fit the data
    # voxs = fit_voxels(x=bold, config=fitting_config)
    res = fit_voxels(x=bold, config=fitting_config)

    # Update the results to match the x0/y0, sigma style used by prfanalyze
    # all_fields = ('index', 'voxel') + fitting_config.fields + ('fit', 'pred',)
    # res = {field_name: voxs[index]
    #        for (index, field_name) in enumerate(all_fields)}

    # Determine the goodness of fit
    r2s = []
    for i in range(len(bold)):
        r2 = pearsonr(res['pred'][i], bold[i])
        if fitting_config.verbose:
            print(f"R2 for voxel {i}: {r2}")
        r2s.append(r2[0])

    # Save results to desired output format and naming
    fitting_config.save_final_results(results=res, shape=bold_im.shape[:-1], affine=bold_im.affine)

    # Also export the prediction and testdata
    fitting_config.save_as_nifti(bold, shape=bold_im.shape, affine=bold_im.affine, filename='testdata.nii.gz')
    # im = nib.Nifti1Image(bold.reshape(bold_im.shape), bold_im.affine)
    # im.to_filename(os.path.join(outdir, 'testdata.nii.gz'))

    fitting_config.save_as_nifti(data=res['pred'], shape=bold_im.shape, affine=bold_im.affine, filename='modelpred.nii.gz')
    # im = nib.Nifti1Image(np.reshape(
    #     res['pred'], bold_im.shape), bold_im.affine)
    # im.to_filename(os.path.join(outdir, 'modelpred.nii.gz'))

    # Store R2
    fitting_config.save_as_nifti(data=r2s, shape=(bold_im.shape[0], 1,1,1), affine=bold_im.affine, filename='r2.nii.gz')
    # im = nib.Nifti1Image(np.reshape(
    #     r2s, (bold_im.shape[0], 1, 1, 1)), bold_im.affine)
    # im.to_filename(os.path.join(outdir, 'r2.nii.gz'))

    # Store config object
    yamls = yaml.YAML()
    # yamls.register_class(FittingConfig)
    with open(os.path.join(outdir, 'options.yml'), 'w') as options_file:
        yamls.dump(opts, options_file)

    # That's it!
    print("PRFPY finished succesfully.")
    sys.exit(0)
