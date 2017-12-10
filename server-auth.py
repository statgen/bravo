from flask import Flask, request, jsonify, abort, url_for, Blueprint, render_template
import requests
from pymongo import MongoClient
import functools
import urllib
import jwt
from datetime import datetime
import hashlib
import os
import argparse


argparser = argparse.ArgumentParser()
argparser.add_argument('--host', default = '0.0.0.0', help = 'the hostname to use to access this server')
argparser.add_argument('--port', type = int, default = 5000, help = 'an integer for the accumulator')


bp = Blueprint('bp', __name__, template_folder = 'templates', static_folder = 'static')

app = Flask(__name__)
app.config.from_object('flask_config.BravoFreeze5GRCh38Config')

BRAVO_AUTH_URL_PREFIX = app.config['BRAVO_AUTH_URL_PREFIX']
BRAVO_AUTH_SECRET = app.config['BRAVO_AUTH_SECRET']
BRAVO_ACCESS_SECRET = app.config['BRAVO_ACCESS_SECRET']

mongo_host = app.config['MONGO']['host']
mongo_port = app.config['MONGO']['port']
mongo_db_name = app.config['MONGO']['name']

GOOGLE_CLIENT_ID = app.config['GOOGLE_LOGIN_CLIENT_ID']
GOOGLE_CLIENT_SECRET = app.config['GOOGLE_LOGIN_CLIENT_SECRET']

GOOGLE_AUTH_API = 'https://accounts.google.com/o/oauth2/v2/auth'
GOOGLE_TOKEN_API = 'https://www.googleapis.com/oauth2/v4/token'
GOOGLE_TOKENINFO_API = 'https://www.googleapis.com/oauth2/v3/tokeninfo'
GOOGLE_AUTH_SCOPE = 'https://www.googleapis.com/auth/userinfo.email'
GOOGLE_ACCESS_TYPE = 'offline'
GOOGLE_RESPONSE_TYPE = 'code'


def setup_auth_tokens_collection(mongo, db_name):
    db = mongo[db_name]
    if 'auth_tokens' not in db.collection_names():
        db.create_collection('auth_tokens')
    db.auth_tokens.create_index('issued_at', expireAfterSeconds = 3 * 60)


mongo = MongoClient(mongo_host, mongo_port, connect = True)
setup_auth_tokens_collection(mongo, mongo_db_name)


def get_db():
   return mongo[mongo_db_name]


def validate_google_access_token(access_token):
    google_response = requests.get(GOOGLE_TOKENINFO_API, params = { 'access_token': access_token })
    if google_response.status_code != 200:
        return (None, None)
    access_token_data = google_response.json()
    expires_in = access_token_data.get('expires_in', None)
    if expires_in is None or int(expires_in) <= 0:
        return (None, None)
    scope = access_token_data.get('scope', None)
    if scope is None or GOOGLE_AUTH_SCOPE not in scope.split():
        return (None, None)
    if access_token_data.get('email_verified', '') != 'true':
        return (None, None)
    return (access_token_data.get('email', None), access_token_data.get('aud', None))


def authorize_user(token_email, token_client_id):
    document = get_db().users.find_one({ 'email': token_email, 'agreed_to_terms': True, 'enabled_api': True }, projection = {'_id': False})
    if not document:
        return False
    if token_client_id != GOOGLE_CLIENT_ID:
        return False
    return True


class UserError(Exception):
    status_code = 400
    def __init__(self, message, status_code = None):
        Exception.__init__(self)
        self.message = message


@bp.errorhandler(UserError)
def handle_user_error(error):
    response = jsonify({ 'error': error.message })
    response.status_code = error.status_code
    return response


@bp.route('/ip', methods = ['GET'])
def ip():
    ip = request.headers.get('X-Forwarded-For', request.remote_addr)
    response = jsonify({ 'ip': ip })
    response.status_code = 200
    return response


@bp.route('/auth', methods = ['GET'])
def auth():
    issued_at = datetime.utcnow()
    ip = request.headers.get('X-Forwarded-For', request.remote_addr)
    auth_token = jwt.encode({'ip': ip, 'iat': issued_at}, BRAVO_AUTH_SECRET, algorithm = 'HS256')
    payload = {
        'client_id': GOOGLE_CLIENT_ID,
        'redirect_uri': url_for('.auth_callback', _external = True, _scheme = 'https'),
        'access_type': GOOGLE_ACCESS_TYPE,
        'response_type': GOOGLE_RESPONSE_TYPE,
        'scope': GOOGLE_AUTH_SCOPE,
        'state': auth_token
    }
    auth_url = '{}?{}'.format(GOOGLE_AUTH_API, urllib.urlencode(payload))
    get_db().auth_tokens.insert({'auth_token': auth_token, 'issued_at': issued_at, 'access_token': None, 'error': None })
    response = jsonify({
        'auth_url': auth_url,
        'auth_token': auth_token
    })
    response.status_code = 200
    return response


@bp.route('/auth/callback', methods = ['GET'])
def auth_callback():
    auth_token = request.args.get('state', None)
    auth_code = request.args.get('code', None)
    try:
        decoded_auth_token = jwt.decode(auth_token, BRAVO_AUTH_SECRET)
    except jwt.InvalidTokenError:
        raise UserError('Bad authorization token.')
    document = get_db().auth_tokens.find_one({ 'auth_token': auth_token, 'access_token': None, 'error': None }, projection = {'_id': False}) 
    if not document:
        raise UserError('Expired authorization token.')
    payload = {
        'client_id': GOOGLE_CLIENT_ID,
        'client_secret': GOOGLE_CLIENT_SECRET,
        'code': auth_code,
        'redirect_uri': url_for('.auth_callback', _external = True, _scheme = 'https'),
        'grant_type': 'authorization_code'
    }
    google_response = requests.post(GOOGLE_TOKEN_API, data = payload)
    google_response_data = google_response.json()
    if google_response.status_code != 200:
        get_db().auth_tokens.update_one({ 'auth_token': auth_token}, {'$set': {'error': google_response_data['error_description']}})
        raise UserError(google_response_data['error_description'])
    email, client_id = validate_google_access_token(google_response_data['access_token'])
    if email is None or client_id is None:
        get_db().auth_tokens.update_one({ 'auth_token': auth_token}, {'$set': {'error': 'Invalid Google access token.'}})
        raise UserError('Invalid Google access token.')
    if not authorize_user(email, client_id):
        get_db().auth_tokens.update_one({ 'auth_token': auth_token}, {'$set': {'error': 'You are not authorized for API access.'}})
        raise UserError('Not authorized')
    access_token = jwt.encode({'email': email, 'ip': decoded_auth_token['ip'], 'iat': datetime.utcnow()}, BRAVO_ACCESS_SECRET, algorithm = 'HS256')
    get_db().auth_tokens.update_one({ 'auth_token': auth_token}, {'$set': {'access_token': access_token}})
    return render_template('auth_ok.html')


@bp.route('/token', methods = ['POST'])
def get_token():
    auth_token = request.form.get('auth_token', None)
    if auth_token is None:
        raise UserError('Bad Request.')
    try:
        decoded_auth_token = jwt.decode(auth_token, BRAVO_AUTH_SECRET)
    except jwt.InvalidTokenError:
        raise UserError('Bad authorization token.')
    ip = request.headers.get('X-Forwarded-For', request.remote_addr)
    if decoded_auth_token['ip'] != ip:
        raise UserError('This authorization token was issued for different IP address.') 
    document = get_db().auth_tokens.find_one({ 'auth_token': auth_token }, projection = {'_id': False})
    if not document:
        raise UserError('Expired authorization token.')
    if document['error'] is not None:
        get_db().auth_tokens.remove({ 'auth_token': auth_token })
        raise UserError(document['error'])
    if document['access_token'] is not None:
        get_db().auth_tokens.remove({ 'auth_token': auth_token })
        response = jsonify({
            'access_token': document['access_token'],
            'token_type': 'Bearer',
            'ip': ip
        })
    else:    
        response = jsonify({
            'access_token': None
        })
    response.status_code = 200
    return response


# We don't check for IP here because we want to allow access token holder be able to revoke access when their IP changed permanently.
@bp.route('/revoke', methods = ['GET'])
def revoke_token():
    access_token = request.args.get('access_token', None)
    if access_token is None:
        raise UserError('Bad Request.')
    try:
        decoded_access_token = jwt.decode(access_token, BRAVO_ACCESS_SECRET)
    except jwt.InvalidTokenError:
        raise UserError('Bad access token.')
    email = decoded_access_token['email']
    issued_at = datetime.utcfromtimestamp(decoded_access_token['iat'])
    document = get_db().users.find_one({ 'email': email }, projection = {'_id': False})
    if not document:
        raise UserError('Bad access token.')
    revoked_at = document.get('access_token_revoked_at', None)
    if revoked_at is not None and revoked_at > issued_at:
        raise UserError('Bad access token.')
    revoked_at = datetime.utcnow()
    get_db().users.update_one({ 'email': email }, {'$set': {'access_token_revoked_at': revoked_at}})
    response = jsonify({
        'revoked': True,
        'revoked_at': revoked_at
    })
    response.status_code = 200
    return response


app.register_blueprint(bp, url_prefix = BRAVO_AUTH_URL_PREFIX)

if __name__ == '__main__':   
    args = argparser.parse_args()
    app.run(host = args.host, port = args.port, threaded = True, use_reloader = True)
