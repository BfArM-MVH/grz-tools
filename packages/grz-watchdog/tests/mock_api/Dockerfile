FROM python:3.12-slim
WORKDIR /app
ENV FLASK_APP=app.py
RUN pip install --no-cache-dir Flask
COPY . .
EXPOSE 5000
CMD ["flask", "run", "--host=0.0.0.0", "--cert=cert.pem", "--key=key.pem"]
