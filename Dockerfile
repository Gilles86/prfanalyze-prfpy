# Use Ubuntu Focal as base
FROM ubuntu:focal-20211006

# Install Node.js, Yarn, and required dependencies
RUN apt-get update && apt-get install -y \
    curl gnupg build-essential git \
    && curl --silent --location https://deb.nodesource.com/setup_16.x | bash - \
    && curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - \
    && echo "deb https://dl.yarnpkg.com/debian/ stable main" | tee /etc/apt/sources.list.d/yarn.list \
    && apt-get update && apt-get install -y nodejs yarn \
    && npm install -g bids-validator \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Conda (Miniforge recommended over deprecated Miniconda2)
RUN curl -L -o ~/miniforge.sh 'https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh' \
    && /bin/bash ~/miniforge.sh -b -p /opt/conda \
    && rm ~/miniforge.sh \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo "conda activate base" >> ~/.bashrc

# Install dependencies for prfpy analysis
RUN . /opt/conda/etc/profile.d/conda.sh && conda create -n prfpy_analysis python=3.9 -y \
    && conda run -n prfpy_analysis pip install --upgrade nilearn nibabel==2.0 h5py wget bids sharedmem pimms \
    setuptools wheel six ruamel.yaml \
    && git clone https://github.com/Gilles86/prfpy.git \
    && conda run -n prfpy_analysis python -m pip install ./prfpy \
    && git clone https://github.com/gallantlab/pycortex.git \
    && conda run -n prfpy_analysis python -m pip install ./pycortex \
    && mkdir -p /root/.config/pycortex/ \
    && touch /root/.config/pycortex/options.cfg

# Set environment variables to ensure conda works in container
ENV PATH="/opt/conda/bin:$PATH"
ENV CONDA_DEFAULT_ENV=prfpy_analysis
ENV CONDA_PREFIX="/opt/conda/envs/prfpy_analysis"
ENV PYTHONPATH=""

# Copy necessary scripts
COPY run.py run_prfpy.py default_config.yml / 
RUN chmod 777 /default_config.yml

COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh


COPY version /version

# Set entrypoint
ENTRYPOINT ["/solve.sh"]

RUN conda run -n prfpy_analysis python -m pip install scikit-learn
RUN conda run -n prfpy_analysis python -m pip install --upgrade --no-cache-dir nilearn
