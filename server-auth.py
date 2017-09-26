from flask import Flask, request, jsonify, abort
import requests
from pymongo import MongoClient
from webargs import fields
from webargs.flaskparser import parser
import functools

app = Flask(__name__)

api_version = 'v1'

port = 7776


GOOGLE_CLIENT_ID = '27789452673-oba39dcqb63aj4q2l0hk9m5dsid366lm.apps.googleusercontent.com'
GOOGLE_CLIENT_SECRET = 'luYXkZfETMgu8Gmyxf-Xwpv-'
GOOGLE_REDIRECT_URI = 'urn:ietf:wg:oauth:2.0:oob'
GOOGLE_TOKEN_API = 'https://www.googleapis.com/oauth2/v4/token'


class UserError(Exception):
    status_code = 400
    def __init__(self, message, status_code = None):
        Exception.__init__(self)
        self.message = message


@app.errorhandler(UserError)
def handle_user_error(error):
    response = jsonify({ 'error': error.message })
    response.status_code = error.status_code
    return response


@app.route('/auth/token', methods = ['POST'])
def get_token():
    auth_code = request.form.get('code', None)
    refresh_token = request.form.get('refresh_token', None)
    if auth_code is not None and refresh_token is not None:
        raise UserError('Bad Request.')
    if auth_code is None and refresh_token is None:
        raise UserError('Bad Request.')
    if auth_code is not None:
        payload = {
            'client_id': GOOGLE_CLIENT_ID,
            'client_secret': GOOGLE_CLIENT_SECRET,
            'code': auth_code,
            'redirect_uri': GOOGLE_REDIRECT_URI,
            'grant_type': 'authorization_code'
        }
        google_response = requests.post(GOOGLE_TOKEN_API, data = payload)
        google_response_data = google_response.json()
        if google_response.status_code != 200:
            raise UserError(google_response_data['error_description'])
        response = jsonify({
            'access_token': google_response_data['access_token'],
            'refresh_token': google_response_data['refresh_token'],
            'token_type': google_response_data['token_type'],
            'expires_in': google_response_data['expires_in']
        })
    elif refresh_token is not None:
        payload = {
            'client_id': GOOGLE_CLIENT_ID,
            'client_secret': GOOGLE_CLIENT_SECRET,
            'refresh_token': refresh_token,
            'grant_type': 'refresh_token'
        }
        google_response = requests.post(GOOGLE_TOKEN_API, data = payload)
        google_response_data = google_response.json()
        if google_response.status_code != 200:
            raise UserError(google_response_data['error_description'])
        response = jsonify({
            'access_token': google_response_data['access_token'],
            'token_type': google_response_data['token_type'],
            'expires_in': google_response_data['expires_in']
        })
    else:
        raise UserError('Invalid request.')
    response.status_code = 200
    return response


if __name__ == '__main__':   
   app.run(host = '0.0.0.0', port = port, debug = True)
