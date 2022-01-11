docker run -i --rm --network=host \
 -v $basedir/BIDS:/bids_dataset \
 -v $basedir:/output \
 -v $basedir/default_config.yml:/config.yml \
 niklasmueller/prfanalyze-prfpy \
 /bids_dataset \
 /config.yml \
 --output_dir /output