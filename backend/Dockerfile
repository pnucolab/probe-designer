FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    gffread \
    samtools \
    wget \
    make \
    pkg-config \
    && cd /tmp \
    && wget -q https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz \
    && tar xzf jellyfish-2.3.1.tar.gz \
    && cd jellyfish-2.3.1 \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd / && rm -rf /tmp/jellyfish-2.3.1* \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY core/ ./core/
COPY backend/ ./backend/
COPY probe_designer.py ./
COPY build_kmer_dbs.py ./
COPY bin/ ./bin/

EXPOSE 8002

CMD ["uvicorn", "backend.main_app:app", "--host", "0.0.0.0", "--port", "8002"]