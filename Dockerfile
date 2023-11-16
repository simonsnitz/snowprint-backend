# Use an official Python image as a parent image with your desired Python version
FROM python:3.10

# Install build dependencies and tools
RUN apt-get update && \
    apt-get install -y build-essential cmake git wget zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Download and compile Diamond Aligner version 2.1.8
RUN wget https://github.com/bbuchfink/diamond/archive/v2.1.8.tar.gz -O diamond.tar.gz && \
    tar -xzvf diamond.tar.gz && \
    mv diamond-2.1.8 diamond && \
    cd diamond && \
    mkdir build && cd build && \
    cmake .. && make && make install && \
    cd ../..

# Set the working directory
WORKDIR /app

# Copy the Python application source code and requirements into the container.
COPY . .

# Install Python dependencies
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    python -m pip install -r requirements.txt

# Expose the port that the application listens on.
EXPOSE 8000

# Run your Python application as the entry point
CMD ["python", "-u", "-B", "main.py"]
