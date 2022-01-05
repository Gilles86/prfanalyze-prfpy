docker run -i --rm --network=host \
 -v $basedir/BIDS:/bids_dataset \
 -v /tank/mueller/projects/analyze-prfpy/default_config.yml:/config.yml \
 niklasmueller/prfanalyze-prfpy \
 /bids_dataset \
 /config.yml
 