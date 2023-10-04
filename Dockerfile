# syntax=docker/dockerfile:1

FROM ubuntu:latest
FROM python:3.10.0

WORKDIR /app

RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    python -m pip install -r requirements.txt

# Copy the source code into the container.
COPY . .

# Expose the port that the application listens on.
EXPOSE 8000

# Run the application.
ENTRYPOINT ["python", "-u", "-B", "main.py"]
