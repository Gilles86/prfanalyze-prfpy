FROM bids/base_fsl

# Install python, nibabel and numpy (nibabel>=2.1 requires python>=3.5, ubuntu trusty has only python 3.4)
RUN apt-get update && \
    apt-get install -y python3

RUN apt-get install python3-pip && \
    pip3 install nibabel==2.0 && \
    # Add all dependencies for analysis here
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
    python setup.py install && \
    cd .. && \
    git clone https://github.com/gallantlab/pycortex.git && \
    cd pycortex && \
    python setup.py install  && \
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
