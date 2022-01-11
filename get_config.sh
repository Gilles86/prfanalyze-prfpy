docker run -i --rm --network=host \
 -v $basedir:/output \
 niklasmueller/prfanalyze-prfpy \
 --get_config /output \