from flask import Flask, request, jsonify, abort, Blueprint
from pymongo import MongoClient
from webargs import fields
from webargs.flaskparser import parser
import functools
import argparse
import jwt
from datetime import datetime
from utils import Xpos
from flask_limiter import Limiter


argparser = argparse.ArgumentParser()
argparser.add_argument('--host', default = '0.0.0.0', help = 'the hostname to use to access this server')
argparser.add_argument('--port', type = int, default = 5000, help = 'an integer for the accumulator')


bp = Blueprint('bp', __name__, template_folder = 'templates', static_folder = 'static')


app = Flask(__name__)
app.config.from_object('flask_config.BravoFreeze5GRCh38Config')

proxy = app.config['PROXY']

BRAVO_API_URL_PREFIX = app.config['BRAVO_API_URL_PREFIX']
BRAVO_ACCESS_SECRET = app.config['BRAVO_ACCESS_SECRET']

mongo_host = app.config['MONGO']['host']
mongo_port = app.config['MONGO']['port']
mongo_db_name = app.config['MONGO']['name']


dataset_name = app.config['DATASET_NAME']
bravo_api_version = app.config['BRAVO_API_VERSION']

pageSize = 10
maxRegion = 100000


mongo = MongoClient(mongo_host, mongo_port, connect = True)


class UserError(Exception):
    status_code = 400
    def __init__(self, message, status_code = None):
        Exception.__init__(self)
        self.message = message


def get_user_ip():
   if proxy:
      x_forwarded_for = request.headers.get('X-Forwarded-For', '').split(',')
      return x_forwarded_for[-1].strip() if len(x_forwarded_for) > 1 else ''
   return request.remote_addr


def get_db():
   return mongo[mongo_db_name]


def validate_access_token(access_token):
   try:
      decoded_access_token = jwt.decode(access_token, BRAVO_ACCESS_SECRET)
   except jwt.InvalidTokenError:
      return (None, None, None)
   return (decoded_access_token.get('email', None), decoded_access_token.get('iat', None), decoded_access_token.get('ip', None))


def authorize_access_token(email, issued_at):
   document = get_db().users.find_one({ 'email': email, 'enabled_api': True, 'agreed_to_terms': True }, projection = {'_id': False})
   if not document:
      return False
   issued_at = datetime.utcfromtimestamp(issued_at)
   revoked_at = document.get('access_token_revoked_at', None)
   if revoked_at is not None and revoked_at > issued_at:
      return False
   return True


def request_is_valid():
   if 'Authorization' not in request.headers:
       return False
   authorization = request.headers.get('Authorization').split()
   if len(authorization) != 2:
      return False
   token_type = authorization[0].lower()
   if token_type != 'bearer':
      return False
   email, issued_at, ip = validate_access_token(authorization[1])
   if email is None or issued_at is None or ip is None:
      return False
   if ip != get_user_ip():
      return False
   return authorize_access_token(email, issued_at)


def require_authorization(func):
   @functools.wraps(func)
   def authorization_wrapper(*args, **kwargs):
      if not request_is_valid():
         raise UserError('not authorized')
      else:
         return func(*args, **kwargs)
   return authorization_wrapper


@parser.error_handler
def handle_parsing_error(error):
   response = jsonify({ 'error': 'invalid query parameters' })
   response.status_code = 400
   abort(response)


@bp.errorhandler(429)
def handle_ratelimit(error):
    response = jsonify({ 'error': 'too many requests' })
    response.status_code = 429
    return response


@bp.errorhandler(UserError)
def handle_user_error(error):
    response = jsonify({ 'error': error.message })
    response.status_code = error.status_code
    return response


@bp.route('/', methods = ['GET'])
@require_authorization
def get_name():
   response = jsonify({
      'dataset': dataset_name,
      'api_version': bravo_api_version
   })
   response.status_code = 200
   return response


@bp.route('/variants', methods = ['GET'])
@require_authorization
def get_variants():
   args = parser.parse({
      'variant_id': fields.Str(required = False),
      'chrom': fields.Str(required = False),
      'start_bp': fields.Int(required = False),
      'end_bp': fields.Int(required = False),
      'position_bp': fields.Int(required = False)
      }, request)
   projection = {'_id': False, 'variant_id': True, 'chrom': True, 'pos': True,  'ref': True, 'alt': True, 'site_quality': True, 'filter': True, 'allele_num': True, 'allele_count': True, 'allele_freq': True}
   if 'variant_id' in args:
      db = get_db()
      try:
         chrom, pos, ref, alt = args['variant_id'].split('-')
         pos = int(pos)
      except ValueError as e:
         raise UserError('Invalid variant name format.')
      if not Xpos.check_chrom(chrom):
         raise UserError('Invalid chromosome name.')
      xpos = Xpos.from_chrom_pos(chrom, pos)
      data = db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection=projection)
   elif all(x in args for x in ['chrom', 'position_bp']):
      db = get_db()
      if not Xpos.check_chrom(args['chrom']):
         raise UserError('Invalid chromosome name.')
      xposition = Xpos.from_chrom_pos(args['chrom'], args['position_bp'])
      data = list(db.variants.find({'xpos': xposition}, projection=projection))
   elif all(x in args for x in ['chrom', 'start_bp', 'end_bp']):
      db = get_db()
      if not Xpos.check_chrom(args['chrom']):
         raise UserError('Invalid chromosome name.')
      xstart = Xpos.from_chrom_pos(args['chrom'], args['start_bp'])
      xend = Xpos.from_chrom_pos(args['chrom'], args['end_bp'])
      if xend - xstart > maxRegion:
         raise UserError('Regions larger than {} base-pairs are not allowed.'.format(maxRegion))
      data = list(db.variants.find({'xpos': {'$lte': xend, '$gte': xstart}}, projection=projection))
   else:
      raise UserError('Invalid Request')
   response = jsonify({
      'data': data
   })
   response.status_code = 200
   return response


app.register_blueprint(bp, url_prefix = BRAVO_API_URL_PREFIX)


limiter = Limiter(app, default_limits = ["180/15 minute"], key_func = get_user_ip)


if __name__ == '__main__':   
   args = argparser.parse_args()
   app.run(host = args.host, port = args.port, threaded = True, use_reloader = True)
