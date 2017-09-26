from flask import Flask, request, jsonify, abort
import requests
from pymongo import MongoClient
from webargs import fields
from webargs.flaskparser import parser
import functools
from rauth import OAuth2Service

app = Flask(__name__)

port = 7776

mongo_host = 'localhost'
mongo_port = 27017
mongo_db_name = 'topmed_freeze5_hg38_testing'

name = 'TOPMed'
version = 'Freeze5'
hg_build = 'GRCh38'
api_version = 'v1'

pageSize = 10
maxRegion = 100000


GOOGLE_CLIENT_ID = '27789452673-oba39dcqb63aj4q2l0hk9m5dsid366lm.apps.googleusercontent.com'
GOOGLE_CLIENT_SECRET = 'luYXkZfETMgu8Gmyxf-Xwpv-'
redirect_uri = 'urn:ietf:wg:oauth:2.0:oob'

google_tokeninfo_api = 'https://www.googleapis.com/oauth2/v3/tokeninfo?access_token={}'
google_authorization_scope = 'https://www.googleapis.com/auth/userinfo.email'
google_token_api = 'https://www.googleapis.com/oauth2/v4/token'


mongo = MongoClient(mongo_host, mongo_port, connect = True)

class Xpos:
    CHROMOSOME_STRINGS = [str(x) for x in range(1, 22+1)] + ['X', 'Y', 'M']
    CHROMOSOME_STRING_TO_NUMBER = {chrom: idx+1 for idx,chrom in enumerate(CHROMOSOME_STRINGS) }
    CHROMOSOME_NUMBER_TO_STRING = {chrom_num: chrom for chrom,chrom_num in CHROMOSOME_STRING_TO_NUMBER.items()}
    @staticmethod
    def from_chrom_pos(chrom, pos):
        if chrom.startswith('chr'): chrom = chrom[3:]
        return Xpos.CHROMOSOME_STRING_TO_NUMBER[chrom] * int(1e9) + pos
    @staticmethod
    def to_chrom_pos(xpos):
        pos = xpos % int(1e9)
        chrom = Xpos.CHROMOSOME_NUMBER_TO_STRING[int(xpos) / int(1e9)]
        return (chrom, pos)
    @staticmethod
    def to_pos(xpos):
        return xpos % int(1e9)
    @staticmethod
    def check_chrom(chrom):
        if chrom.startswith('chr'): chrom = chrom[3:]
        return chrom in Xpos.CHROMOSOME_STRING_TO_NUMBER


class UserError(Exception):
    status_code = 400
    def __init__(self, message, status_code = None):
        Exception.__init__(self)
        self.message = message


def get_db():
   return mongo[mongo_db_name]


# Validates Google OAuth2 access token and returns email and google Client ID
def validate_google_access_token(access_token):
   google_response = requests.get(google_tokeninfo_api.format(access_token))
   if google_response.status_code != 200:
      return (None, None)
   access_token_data = google_response.json()
   if 'aud' not in access_token_data or \
      'expires_in' not in access_token_data or \
      'scope' not in access_token_data or \
      'email' not in access_token_data or \
      'email_verified' not in access_token_data:
      return (None, None)   
   if int(access_token_data['expires_in']) <= 0:
      return (None, None)
   if google_authorization_scope not in access_token_data['scope'].split():
      return (None, None)
   if access_token_data['email_verified'] != 'true':
      return (None, None)
   return (access_token_data['email'], access_token_data['aud'])


# Checks if user with given email exists in our database, has agreed to terms of use, and has registered Google Client ID.
# Returns True if user is authorized for access, otherwise returns False.
def authorize_user(token_email, token_google_client_id):
   document = get_db().users.find_one({ 'email': token_email }, projection = {'_id': False})
   if not document:
      return False
   if not document['agreed_to_terms']:
      return False
   if not document.get('enabled_api', False):
      return False
   user_google_client_id = document.get('google_client_id', None)
   if token_google_client_id != GOOGLE_CLIENT_ID:
      if user_google_client_id is None:
         return false
      if user_google_client_id != google_client_id:
         return False
   return True


def request_is_valid(request):
   if 'Authorization' not in request.headers:
       return False
   authorization = request.headers.get('Authorization').split()
   if len(authorization) != 2:
      return False
   token_type = authorization[0].lower()
   if token_type != 'bearer':
      return False
   email, client_id = validate_google_access_token(authorization[1])
   if not email or not client_id:
      return False
   if not authorize_user(email, client_id):
      return False
   return True


def require_authorization(func):
   @functools.wraps(func)
   def authorization_wrapper(*args, **kwargs):
      if not request_is_valid(request):
         response = jsonify({ 'error': 'not authorized' })
         response.status_code = 400
         return response
      else:
         return func(*args, **kwargs)
   return authorization_wrapper


@parser.error_handler
def handle_parsing_error(error):
   response = jsonify({ 'error': 'invalid query parameters' })
   response.status_code = 400
   abort(response)


@app.errorhandler(UserError)
def handle_user_error(error):
    response = jsonify({ 'error': error.message })
    response.status_code = error.status_code
    return response


@app.route('/', methods = ['GET'])
@require_authorization
def get_name():
   response = jsonify({
      'name': name,
      'version': version,
      'hg_build': hg_build,
      'api_version': api_version
   })
   response.status_code = 200
   return response


@app.route('/variants', methods = ['GET'])
@require_authorization
def get_variants():
   args = parser.parse({
      'chrom': fields.Str(required = True),
      'start': fields.Int(required = True),
      'end': fields.Int(required = True)
      }, request)
   db = get_db()
   if not Xpos.check_chrom(args['chrom']):
      raise UserError('Invalid chromosome name.')
   xstart = Xpos.from_chrom_pos(args['chrom'], args['start'])
   xend = Xpos.from_chrom_pos(args['chrom'], args['end'])
   if xend - xstart > maxRegion:
      raise UserError('Larger than {} bp regions are not allowed.'.format(maxRegion))
   variants = list(db.variants.find({'xpos': {'$lte': xend, '$gte': xstart}}, projection={'_id': False}))
   response = jsonify({
      'data': variants
   })
   response.status_code = 200
   return response


if __name__ == '__main__':   
   app.run(host = '0.0.0.0', port = port, debug = True)
