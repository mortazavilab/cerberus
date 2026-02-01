# ---------- build stage ----------
FROM continuumio/miniconda3:latest AS builder

RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    ca-certificates \
    openssl \
    pkg-config \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /wheels

RUN printf "\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
channel_priority: strict\n\
ssl_verify: false\n" > /opt/conda/.condarc

ENV CONDA_SSL_VERIFY=no

RUN conda create -y -n app_env python=3.8
RUN conda run -n app_env conda install -y -c conda-forge \
  pytables=3.8.0 \
  h5py=3.11.0 \
  hdf5
RUN conda run -n app_env pip install --no-cache-dir --upgrade pip wheel

RUN conda run -n app_env pip wheel \
    --wheel-dir=/wheels \
    streamlit==1.40.1 \
    numpy==1.24.4 \
    pandas==1.5.3 \
    anndata==0.9.2 \
    scanpy==1.9.8 \
    scipy==1.10.1 \
    pyranges==0.1.4 \
    plotly \
    git+https://github.com/fairliereese/python-ternary \
    git+https://github.com/mortazavilab/swan_vis.git \
    git+https://github.com/mortazavilab/cerberus.git

# ------------ runtime stage---------------
FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y \
    ca-certificates \
    openssl \
    libblas3 \
    liblapack3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /opt/conda/envs/app_env /opt/conda/envs/app_env
COPY --from=builder /wheels /wheels

RUN printf "\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
channel_priority: strict\n\
ssl_verify: false\n" > /opt/conda/.condarc

ENV CONDA_SSL_VERIFY=no

RUN conda run -n app_env pip install --no-cache-dir /wheels/* \
    && rm -rf /wheels \
    && rm -rf /opt/conda/pkgs/*

WORKDIR /app
COPY *.py /app/
COPY *.png /app/
COPY tabs/*.py /app/tabs/

SHELL ["conda", "run", "-n", "app_env", "/bin/bash", "-c"]

EXPOSE 8501
CMD ["conda", "run", "--no-capture-output", "-n", "app_env", "streamlit", "run", "/app/main.py", "--server.port=8501", "--server.address=0.0.0.0"]
