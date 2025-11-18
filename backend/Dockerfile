FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    gffread \
    samtools \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY core/ ./core/
COPY backend/ ./backend/
COPY simple_pipeline.py ./
COPY bin/ ./bin/

EXPOSE 8002

CMD ["uvicorn", "backend.main_app:app", "--host", "0.0.0.0", "--port", "8002"]