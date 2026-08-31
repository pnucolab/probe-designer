# GPU-enabled image for the probe pipeline.
#
# The aligner and k-mer counting run on the GPU via CuPy, which JIT-compiles
# CUDA kernels — so the image needs the CUDA toolkit (NVRTC + headers), not just
# the runtime. Jellyfish and RazerS3 remain for the ungapped/edge fallbacks.
#
# The container only accelerates when GPUs are attached: grant them in
# docker-compose.yml (deploy.resources.reservations.devices) and install
# nvidia-container-toolkit on the host. Without a GPU the pipeline still runs
# (RazerS3 + Jellyfish) but without acceleration.
FROM nvidia/cuda:12.9.1-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    CUDA_PATH=/usr/local/cuda \
    PATH=/usr/local/cuda/bin:/usr/local/bin:${PATH}

WORKDIR /app

# Python 3.11 (deadsnakes), build tools, samtools/gffread, and Jellyfish (from source).
RUN apt-get update && apt-get install -y --no-install-recommends \
        software-properties-common ca-certificates curl \
    && add-apt-repository -y universe \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get update && apt-get install -y --no-install-recommends \
        python3.11 python3.11-dev python3.11-distutils \
        gcc g++ git samtools gffread wget make pkg-config \
    && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11 \
    && ln -sf /usr/bin/python3.11 /usr/local/bin/python \
    && ln -sf /usr/bin/python3.11 /usr/local/bin/python3 \
    && cd /tmp \
    && wget -q https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz \
    && tar xzf jellyfish-2.3.1.tar.gz && cd jellyfish-2.3.1 \
    && ./configure --prefix=/usr/local && make -j"$(nproc)" && make install && ldconfig \
    && cd / && rm -rf /tmp/jellyfish-2.3.1* /var/lib/apt/lists/*

COPY requirements.txt .
RUN python -m pip install --no-cache-dir -r requirements.txt \
    && python -m pip install --no-cache-dir cupy-cuda12x==13.5.1

COPY core/ ./core/
COPY backend/ ./backend/
COPY config/ ./config/
COPY probe_designer.py ./
COPY build_kmer_dbs.py ./

# NOTE: `COPY bin/ ./bin/` was removed. The pipeline is GPU-only and no longer
# shells out to RazerS3, and gffread is installed from apt above, so nothing in
# the image ever read /app/bin. The two prebuilt binaries it held are not
# tracked in this repository.

EXPOSE 8002

CMD ["uvicorn", "backend.main_app:app", "--host", "0.0.0.0", "--port", "8002"]
