from flask import Flask, request, jsonify, abort, Blueprint
from pymongo import MongoClient, ASCENDING, DESCENDING
from webargs import fields, ValidationError
from webargs.flaskparser import parser
import functools
import argparse
import jwt
import string
import bson
from datetime import datetime
from utils import Xpos
from flask_limiter import Limiter
import os

argparser = argparse.ArgumentParser()
argparser.add_argument('--host', default = '0.0.0.0', help = 'The hostname to use to access this server.')
argparser.add_argument('--port', type = int, default = 5000, help = 'An integer for the port number.')

bp = Blueprint('bp', __name__, template_folder = 'templates', static_folder = 'static')

app = Flask(__name__, instance_relative_config = True)

# Load default config
app.config.from_object('config.default')
# Load instance configuration if exists
app.config.from_pyfile('config.py', silent = True)
# Load configuration file specified in BRAVO_CONFIG_FILE environment variable if exists
app.config.from_envvar('BRAVO_CONFIG_FILE', silent = True)


proxy = app.config['PROXY']

URL_PREFIX = app.config['URL_PREFIX']
API_URL_PREFIX = app.config['API_URL_PREFIX']
BRAVO_ACCESS_SECRET = app.config['BRAVO_ACCESS_SECRET']

mongo_host = app.config['MONGO']['host']
mongo_port = app.config['MONGO']['port']
mongo_db_name = app.config['MONGO']['name']

api_dataset_name = app.config['API_DATASET_NAME']
api_collection_name = app.config['API_COLLECTION_NAME']
api_version = app.config['API_VERSION']

pageSize = app.config['API_PAGE_SIZE']
maxRegion = app.config['API_MAX_REGION']

projection = {'_id': True, 'xpos': True, 'variant_id': True, 'chrom': True, 'pos': True,  'ref': True, 'alt': True, 'site_quality': True, 'filter': True, 'allele_num': True, 'allele_count': True, 'allele_freq': True, 'rsids': True}
allowed_sort_keys = {'pos': long, 'allele_count': int, 'allele_freq': float, 'allele_num': int, 'site_quality': float, 'filter': str, 'variant-id': str}
allowed_filter_keys = {'allele_count', 'allele_freq', 'allele_num', 'site_quality', 'filter'}

vcf_header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
vcf_meta = [
         '##fileformat=VCFv4.2',
         '##FILTER=<ID=PASS,Description="All filters passed">',
         '##FILTER=<ID=CEN,Description="Variant located in centromeric region with inferred sequences">',
         '##FILTER=<ID=SVM,Description="Variant failed SVM filter">',
         '##FILTER=<ID=DISC,Description="Mendelian or duplicate genotype discordance is high (3/5% or more)">',
         '##FILTER=<ID=CHRXHET,Description="Excess heterozygosity in chrX in males">',
         '##FILTER=<ID=EXHET,Description="Excess heterozygosity with HWE p-value < 1e-6">',
         '##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of Alleles in Samples with Coverage">',
         '##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate Allele Counts in Samples with Coverage">',
         '##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate Allele Frequencies">'
      ]

mongo = MongoClient(mongo_host, mongo_port, connect = True)


class UserError(Exception):
    status_code = 400
    def __init__(self, message, status_code = None):
        Exception.__init__(self)
        self.message = message


def get_user_ip():
   if proxy:
      x_forwarded_for = request.headers.get('X-Forwarded-For', '').split(',')
      return x_forwarded_for[-1].strip() if len(x_forwarded_for) > 0 else ''
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
   if app.config['API_GOOGLE_AUTH']:
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
   elif app.config['API_IP_WHITELIST']:
      if get_user_ip() in app.config['API_IP_WHITELIST']:
         return True
   else:
      return True


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
      'dataset': api_dataset_name,
      'api_version': api_version
   })
   response.status_code = 200
   return response


@bp.route('/variant', methods = ['GET'])
@require_authorization
def get_variant():
   args = parser.parse({
      'variant_id': fields.Str(required = False, validate = lambda x: len(x) > 0),
      'chrom': fields.Str(required = False, validate = lambda x: len(x) > 0),
      'pos': fields.Int(required = False, validate = lambda x: x >= 0),
      'vcf': fields.Bool(required = False, missing = False)
      }, request)

   if 'variant_id' in args:
      if args['variant_id'].startswith('rs'):
          mongo_filter = { 'rsids' : args['variant_id'] }
      else:
          try:
              chrom, pos, ref, alt = args['variant_id'].split('-')
              pos = int(pos)
          except ValueError as e:
              raise UserError('Invalid variant name format.')
          if not Xpos.check_chrom(chrom):
              raise UserError('Invalid chromosome name.')
          xpos = Xpos.from_chrom_pos(chrom, pos)
          mongo_filter = { 'xpos': xpos, 'ref': ref, 'alt': alt }
   elif all(x in args for x in ['chrom', 'pos']):
      if not Xpos.check_chrom(args['chrom']):
         raise UserError('Invalid chromosome name.')
      xpos = Xpos.from_chrom_pos(args['chrom'], args['pos'])
      mongo_filter = { 'xpos': xpos }
   else:
      raise UserError('Invalid query parameters.')

   data = []
   response = { 'next': None }
   db = get_db()
   collection = db[api_collection_name]
   cursor = collection.find(mongo_filter, projection=projection)
   if not args['vcf']:
      response['format'] = 'json'
      for r in cursor:
         last_object_id = r.pop('_id')
         r.pop('xpos', None)
         data.append(r)
         last_variant = r
   else:
      response['format'] = 'vcf'
      response['header'] = vcf_header
      response['meta'] = vcf_meta
      for r in cursor:
         last_object_id = r.pop('_id')
         r.pop('xpos', None)
         data.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\tAN={};AC={};AF={}'.format(
            r['chrom'], r['pos'], ';'.join(r['rsids']) if r['rsids'] else '.', r['ref'], r['alt'], r['site_quality'], r['filter'],
            r['allele_num'], r['allele_count'], r['allele_freq']))
         last_variant = r
   response['data'] = data

   response = jsonify(response)
   response.status_code = 200
   return response



def deserialize_query_sort(value):
   query_sort = list()
   for key_direction in (x.strip().split(':') for x in value.strip().split(',')):
      if len(key_direction) == 0 or len(key_direction) > 2:
         raise ValidationError('invalid sort syntax')
      key = key_direction[0].lower()
      if key not in allowed_sort_keys:
         raise ValidationError('sort is not supported for specified key(s)')
      direction = key_direction[1].lower() if len(key_direction) > 1 else ''
      if direction == '' or direction == 'asc' or direction == '1':
         direction = ASCENDING
      elif direction == 'desc' or direction == '-1':
         direction = DESCENDING
      else:
         raise ValidationError('unkown sort direction')
      query_sort.append((key, direction))
   return query_sort


def deserialize_query_filter(value, value_type):
   value = value.strip()
   if len(value) == 0:
      raise ValidationError('empty value')
   if len(value) > 1 and ((value[0] == '"' and value[-1] == '"') or (value[0] == '\'' and value[-1] == '\'')):
      operator = '$eq'
      value = value[1:-1]
      return { operator: value }
   elements = value.split(':')
   if len(elements) == 1:
      operator = '$eq'
      value = elements[0]
   elif len(elements) == 2:
      operator = '${}'.format(elements[0].lower())
      value = elements[1]
      if operator not in ['$eq', '$ne', '$gt', '$lt', '$gte', '$lte']:
         raise ValidationError('unknown operator')
   else:
      raise ValidationError('too many values')
   try:
      if value_type == int:
         value = int(value)
      elif value_type == long:
         value = long(value)
      elif value_type == float:
         value = float(value)
      elif value_type == str:
         if len(value) > 1 and ((value[0] == '"' and value[-1] == '"') or (value[0] == '\'' and value[-1] == '\'')):
            value = value[1:-1]
      else:
         raise ValidationError('unsupported value type')
   except:
      raise ValidationError('invalid value type')
   return (operator, value)


def deserialize_query_last(value):
   elements = value.strip().split(':')
   if len(elements) == 0:
      raise ValidationError('empty value')
   objectid = elements[0]
   if len(objectid) != 24 or any(c not in string.hexdigits for c in objectid):
      raise ValidationError('invalid value')
   return elements


def validate_query(value):
   # eliminate duplicated entries in filters and sort
   for key in allowed_filter_keys:
      key_filters = value.get(key, None)
      if key_filters is not None:
         value[key] = dict((x for x in key_filters))
   if 'sort' in value:
      value['sort'] = dict(value['sort']).items()
   # check if last element is consistent with sort fields and cast to the corresponding types
   if 'last' in value and len(value['last']) > 1:
      if 'sort' not in value or len(value['last']) - 1 != len(value['sort']):
         return False
      try:
         for i, key in enumerate((x for x, y in value['sort']), 1):
            key_type = allowed_sort_keys[key]
            if key_type == int:
               value['last'][i] = int(value['last'][i])
            elif key_type == long:
               value['last'][i] = long(value['last'][i])
            elif key_type == float:
               value['last'][i] = float(value['last'][i])
      except:
         return False
   return True


def build_region_query(args, xstart, xend):
   # prepare sort conditions in mongo format
   if 'sort' not in args or len(args['sort']) == 0:
      mongo_sort = [(u'xpos', ASCENDING)] # xpos sorted by default if nothing else is specified
   else:
      mongo_sort = [(u'xpos', x[1]) if x[0] == 'pos' else x for x in args['sort']] # if user was sorting by 'pos', then replace 'pos' to 'xpos'

   mongo_filter = {
      'xpos' : {'$lte': xend, '$gte': xstart}
   }

   for filter_key in allowed_filter_keys:
      filter_values = args.get(filter_key, None)
      if filter_values is not None:
         mongo_filter[filter_key] = filter_values

   if 'last' in args:
      objectid = args['last'][0]

      if len(mongo_sort) > 0:
         mongo_filter['$or'] = list()

         key, sort_direction = mongo_sort[0]
         last_value = args['last'][1]
         if key not in mongo_filter:
            mongo_filter[key] = {}
         if sort_direction == ASCENDING:
            mongo_filter[key].update({ '$gte': last_value })
            mongo_filter['$or'].append({key: {'$ne': last_value}})
         else:
            mongo_filter[key].update({ '$lte': last_value })
            mongo_filter['$or'].append({key: {'$ne': last_value}})

         if len(mongo_sort) > 1:
            mongo_filter['$or'].append({ '$and':[ { '$or': [ {'_id': {'$gt': bson.objectid.ObjectId(objectid)}}  ] } ]})
            for i in xrange(len(mongo_sort) - 1, 0, -1):
               key, sort_direction = mongo_sort[i]
               last_value = args['last'][i + 1]
               for filters in mongo_filter['$or'][-1]['$and']:
                  filters['$or'].append({key: {'$ne': last_value}})
               mongo_filter['$or'][-1]['$and'].append({ '$or': [ {key: { '$gte' if sort_direction == ASCENDING else '$lte': last_value}} ]})
         else:
            mongo_filter['$or'].append({'_id': {'$gt': bson.objectid.ObjectId(objectid)}})
      else:
         mongo_filter['_id'] = {'$gt': bson.objectid.ObjectId(objectid)}
   return mongo_filter, mongo_sort


def build_link_next(args, last_object_id, last_variant, mongo_sort):
   if last_object_id is None or last_variant is None:
      return None
   link_next = request.base_url + '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems() if arg != 'last'))
   if 'sort' not in args or len(args['sort']) == 0:
      if len(mongo_sort) > 0:
         link_next += '&sort=' +  ','.join(('{}:{}'.format('pos' if key == 'xpos' else key, 'asc' if sort_direction == ASCENDING else 'desc') for key, sort_direction in mongo_sort))
   link_next = '{}&last={}'.format(link_next, last_object_id)
   for key, sort_direction in mongo_sort:
      link_next += ':{}'.format(last_variant[key] if key != 'xpos' else Xpos.from_chrom_pos(last_variant['chrom'], last_variant['pos']))
   return link_next


@bp.route('/region', methods = ['GET'])
@require_authorization
def get_region():
   arguments = {
       'chrom': fields.Str(required = True, validate = lambda x: len(x) > 0),
       'start': fields.Int(required = True, validate = lambda x: x >= 0),
       'end': fields.Int(required = True, validate = lambda x: x > 0),
       'allele_count': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, int))),
       'allele_freq': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, float))),
       'allele_num': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, int))),
       'site_quality': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, float))),
       'filter': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, str))),
       'sort': fields.Function(deserialize = deserialize_query_sort),
       'vcf': fields.Bool(required = False, missing = False),
       'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = pageSize),
       'last': fields.Function(deserialize = deserialize_query_last)
   }

   args = parser.parse(arguments, validate = validate_query)

   #print args

   if args['start'] >= args['end']:
      raise UserError('Start position must be less than end position.')

   # check if start less than end
   if args['end'] - args['start'] > maxRegion:
      raise UserError('Regions larger than {} base-pairs are not allowed.'.format(maxRegion))

   # transform start and end to xpos
   if not Xpos.check_chrom(args['chrom']):
      raise UserError('Invalid chromosome name.')
   xstart = Xpos.from_chrom_pos(args['chrom'], args['start'])
   xend = Xpos.from_chrom_pos(args['chrom'], args['end'])

   mongo_filter, mongo_sort = build_region_query(args, xstart, xend)
   # print mongo_filter
   # print mongo_sort

   data = [];
   response = {}
   last_variant = None
   last_object_id = None
   db = get_db()
   collection = db[api_collection_name]
   cursor = collection.find(mongo_filter, projection = projection).sort(mongo_sort).limit(args['limit']) 
   if not args['vcf']:
      response['format'] = 'json'
      for r in cursor:
         last_object_id = r.pop('_id')
         r.pop('xpos', None)
         data.append(r)
         last_variant = r
   else:
      response['format'] = 'vcf'
      response['header'] = vcf_header
      response['meta'] = vcf_meta
      for r in cursor:
         last_object_id = r.pop('_id')
         r.pop('xpos', None)
         data.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\tAN={};AC={};AF={}'.format(
            r['chrom'], r['pos'], ';'.join(r['rsids']) if r['rsids'] else '.', r['ref'], r['alt'], r['site_quality'], r['filter'],
            r['allele_num'], r['allele_count'], r['allele_freq']))
         last_variant = r

   response['data'] = data
   response['next'] = build_link_next(args, last_object_id, last_variant, mongo_sort) if len(data) == args['limit'] else None
   response = jsonify(response)
   response.status_code = 200
   return response


@bp.route('/gene', methods = ['GET'])
@require_authorization
def get_gene():
   arguments = {
       'name': fields.Str(required = True, validate = lambda x: len(x) > 0),
       'allele_count': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, int))),
       'allele_freq': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, float))),
       'allele_num': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, int))),
       'site_quality': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, float))),
       'filter': fields.List(fields.Function(deserialize = lambda x: deserialize_query_filter(x, str))),
       'sort': fields.Function(deserialize = deserialize_query_sort),
       'vcf': fields.Bool(required = False, missing = False),
       'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = pageSize),
       'last': fields.Function(deserialize = deserialize_query_last)
   }

   args = parser.parse(arguments, validate = validate_query)

   db = get_db()
   gene = db.genes.find_one({'$or': [
         {'gene_id': args['name']}, {'gene_name': args['name']}, {'other_names': args['name']}
      ]}, projection={'_id': False})
   if not gene:
      raise UserError('Gene with name or identifier equal to {} was not found.'.format(args['name']))

   response = {
      'gene': {
         'name': gene['gene_name'],
         'id': gene['gene_id'],
         'chrom': gene['chrom'],
         'start': gene['start'],
         'stop': gene['stop'],
         'strand': gene['strand']
      }
   }

   mongo_filter, mongo_sort = build_region_query(args, gene['xstart'], gene['xstop'])

   data = [];
   last_variant = None
   last_object_id = None
   db = get_db()
   collection = db[api_collection_name]
   cursor = collection.find(mongo_filter, projection = projection).sort(mongo_sort).limit(args['limit'])
   if not args['vcf']:
      response['format'] = 'json'
      for r in cursor:
         last_object_id = r.pop('_id')
         r.pop('xpos', None)
         data.append(r)
         last_variant = r
   else:
      response['format'] = 'vcf'
      response['header'] = vcf_header
      response['meta'] = vcf_meta
      for r in cursor:
         last_object_id = r.pop('_id')
         r.pop('xpos', None)
         data.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\tAN={};AC={};AF={}'.format(
            r['chrom'], r['pos'], ';'.join(r['rsids']) if r['rsids'] else '.', r['ref'], r['alt'], r['site_quality'], r['filter'],
            r['allele_num'], r['allele_count'], r['allele_freq']))
         last_variant = r

   response['data'] = data
   response['next'] = build_link_next(args, last_object_id, last_variant, mongo_sort) if len(data) == args['limit'] else None
   response = jsonify(response)
   response.status_code = 200
   return response


limiter = Limiter(app, default_limits = app.config['API_REQUESTS_RATE_LIMIT'], key_func = get_user_ip)


app.register_blueprint(bp, url_prefix = URL_PREFIX + API_URL_PREFIX)


if __name__ == '__main__':   
   args = argparser.parse_args()
   app.run(host = args.host, port = args.port, threaded = True, use_reloader = True)
