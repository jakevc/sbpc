# Build stage
FROM rust:slim-bullseye as builder

# Install dependencies
RUN apt-get update && apt-get install -y \
    libhts-dev \
    libcurl4-gnutls-dev \
    libdeflate-dev \
    liblzma-dev \
    zlib1g-dev \
    cmake \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Create a new empty project
WORKDIR /app
COPY . .

# Build the project
RUN cargo build --release

# Runtime stage
FROM debian:bullseye-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libhts3 \
    libcurl4-gnutls-dev \
    libdeflate0 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /app/target/release/sbpc /usr/local/bin/sbpc

# Set entrypoint
ENTRYPOINT ["sbpc"]
CMD ["--help"]
