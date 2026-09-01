from flask import Flask, jsonify

app = Flask("bfarm-mock-api")


@app.route("/ping")
def ping():
    return jsonify({"ping": "pong"})


@app.route("/token", methods=["POST", "GET"])
def get_token():
    return jsonify({"access_token": "mock-token-12345", "expires_in": 3600}), 200


@app.route("/api/v1/upload", methods=["POST", "PUT"])
def upload_pruefbericht():
    return jsonify({"status": "success", "message": "Pruefbericht received by mock API."}), 200


@app.route("/path:path", methods=["POST", "GET", "PUT", "DELETE"])
def catch_all(path):
    return jsonify({"message": f"Endpoint '/{path}' mocked."}), 200
