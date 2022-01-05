FROM ubuntu:focal-20211006

# Install python, nibabel and numpy (nibabel>=2.1 requires python>=3.5, ubuntu trusty has only python 3.4)
RUN apt-get update  \
    && apt-get install -y python3 \
# Install Node.js, Yarn and required dependencies
    && apt-get install -y curl gnupg build-essential \
    && curl --silent --location https://deb.nodesource.com/setup_16.x | bash - \
    && curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - \
    && echo "deb https://dl.yarnpkg.com/debian/ stable main" | tee /etc/apt/sources.list.d/yarn.list \
    && apt-get remove -y --purge cmdtest \
    && apt-get update \
    && apt-get install -y nodejs yarn \
    # remove useless files from the current layer
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /var/lib/apt/lists.d/* \
    && apt-get autoremove \
    && apt-get clean \
    && apt-get autoclean  

# RUN npm install --global npm@^7
RUN npm install -g bids-validator

RUN apt-get update && \
    apt-get install -y git python3-pip && \
    # Add all dependencies for analysis here
    pip3 install nibabel==2.0 && \
    pip3 install --upgrade numpy && \
    pip3 install --upgrade nilearn && \
    pip3 install --upgrade h5py && \
    pip3 install --upgrade wget && \
    pip3 install --upgrade bids && \
    pip3 install --upgrade sharedmem && \
    pip3 install --upgrade pimms && \
    pip3 install -U setuptools wheel && \
    git clone https://github.com/VU-Cog-Sci/prfpy.git && \
    cd prfpy && \
    python3 setup.py install && \
    cd .. && \
    git clone https://github.com/gallantlab/pycortex.git && \
    cd pycortex && \
    python3 setup.py install  && \
    cd .. && \
    mkdir -p /root/.config/pycortex/ && \
    touch /root/.config/pycortex/options.cfg && \
    # 
    apt-get remove -y python3-pip && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



# RUN . /opt/conda/etc/profile.d/conda.sh && \
#     conda create -n prfpy_analysis --channel intel intelpython3_full && \
#     conda activate prfpy_analysis && \
#     pip install --upgrade nilearn && \
#     pip install --upgrade nibabel && \
#     pip install --upgrade h5py && \
#     pip install --upgrade wget && \
#     pip install --upgrade bids && \
#     pip install --upgrade sharedmem && \
#     pip install --upgrade pimms && \
#     pip install -U setuptools wheel && \
#     git clone https://github.com/VU-Cog-Sci/prfpy.git && \
#     cd prfpy && \
#     python setup.py install && \
#     cd .. && \
#     git clone https://github.com/gallantlab/pycortex.git && \
#     cd pycortex && \
#     python setup.py install  && \
#     cd .. && \
#     mkdir -p /root/.config/pycortex/ && \
#     touch /root/.config/pycortex/options.cfg 

ENV PYTHONPATH=""

COPY run.py /run.py
COPY run_prfpy.py /run_prfpy.py

COPY version /version

ENTRYPOINT ["/run.py"]
