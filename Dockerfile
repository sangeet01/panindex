# FractalIndex (FRX) - Docker Image
#
# Build:
#   docker build -t frx .
#
# Run CLI (default):
#   docker run frx frx --help
#   docker run -v $(pwd)/data:/app/data frx frx stats --index /app/data/index.db
#
# Run API server:
#   docker run -p 5000:5000 -e FRX_API_TOKEN=mysecret frx python viz/server.py
#
# Or with docker compose:
#   docker compose up frx-api

FROM python:3.10-slim

LABEL maintainer="FractalIndex Project"
LABEL description="Content-addressed pangenome indexing with alignment-free query"

WORKDIR /app

# Install system dependencies (none needed - pure Python)
RUN apt-get update && apt-get install -y --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies first (better layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source
COPY src/ ./src/
COPY viz/ ./viz/
COPY pyproject.toml .

# Install the package (editable for simplicity in dev; use --no-deps since
# requirements.txt already installed everything)
RUN pip install --no-cache-dir -e . --no-deps

# Create non-root user for security
RUN useradd -m -u 1000 frxuser && chown -R frxuser:frxuser /app
USER frxuser

# Create a data directory for mounted volumes
RUN mkdir -p /app/data

EXPOSE 5000

# Environment variables
ENV FRX_API_TOKEN=""
ENV FRX_RATE_LIMIT="100"
ENV PYTHONUNBUFFERED=1

# Health check using the public /health endpoint
HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 \
  CMD python -c \
      "import urllib.request, sys; \
       r = urllib.request.urlopen('http://localhost:5000/health', timeout=3); \
       sys.exit(0 if r.status == 200 else 1)"

# Default entrypoint is the frx CLI
# Override CMD to run the API server:
#   docker run frx python viz/server.py
ENTRYPOINT ["python", "-m", "src.cli"]
CMD ["--help"]
