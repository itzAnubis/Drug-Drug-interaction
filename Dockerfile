# Use Python 3.11 (best for transformers/RDKit balance)
FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libsm6 \
    libxext6 \
    g++ \
    build-essential \
    cmake && \
    rm -rf /var/lib/apt/lists/*

# Set up virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies with version pins
COPY requirements.txt .
RUN pip install --no-cache-dir \
    numpy==1.26.4 && \
    pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir sentencepiece

# Copy application
COPY . .

# Run FastAPI
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]