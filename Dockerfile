FROM ubuntu:focal-20211006

# Install Node.js, Yarn and required dependencies
RUN apt-get update  \
    && apt-get install -y curl gnupg build-essential \
    && curl --silent --location https://deb.nodesource.com/setup_16.x | bash - \
    && curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - \
    && echo "deb https://dl.yarnpkg.com/debian/ stable main" | tee /etc/apt/sources.list.d/yarn.list \
    && apt-get remove -y --purge cmdtest \
    && apt-get update \
    && apt-get install -y nodejs yarn \
    && apt-get install -y git \
    # remove useless files from the current layer
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /var/lib/apt/lists.d/* \
    && apt-get autoremove \
    && apt-get clean \
    && apt-get autoclean  \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN npm install -g bids-validator

# Install conda
RUN curl -L -o ~/miniconda.sh 'https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh'
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda \
 && rm ~/miniconda.sh \
 && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
 && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
 && echo "conda activate base" >> ~/.bashrc

# Install dependecies for prfpy analysis
RUN . /opt/conda/etc/profile.d/conda.sh && \
    conda create -n prfpy_analysis --channel intel intelpython3_full && \
    conda activate prfpy_analysis && \
    pip install --upgrade nilearn && \
    pip install nibabel==2.0 && \
    pip install --upgrade h5py && \
    pip install --upgrade wget && \
    pip install --upgrade bids && \
    pip install --upgrade sharedmem && \
    pip install --upgrade pimms && \
    pip install -U setuptools wheel && \
    pip install --upgrade -U six && \
    pip install --upgrade -U ruamel.yaml && \
    git clone https://github.com/VU-Cog-Sci/prfpy.git && \
    cd prfpy && \
    python setup.py install && \
    cd .. && \
    git clone https://github.com/gallantlab/pycortex.git && \
    cd pycortex && \
    python setup.py install  && \
    cd .. && \
    mkdir -p /root/.config/pycortex/ && \
    touch /root/.config/pycortex/options.cfg 

ENV PYTHONPATH=""

COPY run.py /run.py
COPY run_prfpy.py /run_prfpy.py
COPY default_config.yml /default_config.yml

COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh

COPY version /version

ENTRYPOINT ["/solve.sh"]
