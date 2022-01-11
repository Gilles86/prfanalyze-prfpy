## PRF-Analyze PRFPY
This app is part of a validation framework for fMRI analysis tools. It builds upon the BIDS file structure and can be used to verify the integrity in terms of reproduction of ground truth of e.g. prf-models.

This app is an alternative implementation to ![this](https://github.com/vistalab/prfmodel) validation framework, i.e. the prfanalyze part is taking over here. Hence, the prfsynth dockerimage should be used to synthesize data in order to know the ground truth. After using this app, the prfreport dockerimage can be used to create meaningful and insightful visual reports about the goodness of fit achieved by this and/or other analysis tools.

### Description
This app will apply a user-specified analysis to a given BIDS dataset.

<!-- ### Documentation
Provide a link to the documentation of your pipeline. -->

<!-- ### How to report errors
Provide instructions for users on how to get help and report errors. -->

<!-- ### Acknowledgments
Describe how would you would like users to acknowledge use of your App in their papers (citation, a paragraph that can be copy pasted, etc.) -->

### Usage
This App has the following command line arguments:

		usage: run.py [-h]
		              [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
		              bids_dir config_file

		prfanalyze-prfpy entry point script.

		positional arguments:
		  bids_dir              The directory with the input dataset formatted
		                        according to the BIDS standard.
		  config_file			The path to a yaml file containing all necessary information to perform the 							  prf-fitting

		optional arguments:
		  -h, --help            show this help message and exit
		  --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
		                        The label(s) of the participant(s) that should be
		                        analyzed. The label corresponds to
		                        sub-<participant_label> from the BIDS spec (so it does
		                        not include "sub-"). If this parameter is not provided
		                        all subjects will be analyzed. Multiple participants
		                        can be specified with a space separated list.

<!-- output_dir            The directory where the output files should be stored.
					If you are running a group level analysis, this folder
					should be prepopulated with the results of
					the participant level analysis. -->
					
To run it in (for one participant):

	docker run -i --rm \
 		-v $basedir/BIDS:/bids_dataset \
 		-v analyze-prfpy/default_config.yml:/config.yml \
 		niklasmueller/prfanalyze-prfpy \
 		/bids_dataset \
 		/config.yml --participant_label 01
 
<!-- 
### Special considerations
Describe whether your app has any special requirements. For example:

- Multiple map reduce steps (participant, group, participant2, group2 etc.)
- Unusual memory requirements
- etc. -->
