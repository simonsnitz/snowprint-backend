# Use an official Python image as a parent image with your desired Python version
FROM python:3.10

# Install build dependencies and tools
RUN apt-get update && \
    apt-get install -y wget build-essential \
    && wget https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz \
    && tar xzf diamond-linux64.tar.gz \
    && mv diamond /usr/bin/ \
    && rm diamond-linux64.tar.gz

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
