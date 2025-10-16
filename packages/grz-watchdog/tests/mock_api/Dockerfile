FROM docker.io/python:3.14-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && \
    rm -rf /var/lib/apt/lists/*
RUN pip install Flask==3.1.2

WORKDIR /app
COPY ./*.py .

EXPOSE 5000

CMD ["flask", "--app", "app", "run", "--host=0.0.0.0", "--port=5000", "--cert=cert.pem", "--key=key.pem"]
