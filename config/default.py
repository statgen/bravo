DEBUG = False
TESTING = False
PROXY = False # True if app is proxied by Apache or similar. 

BROWSER_NAME = 'Bravo'
DATASET_NAME = 'Example Dataset'
SHOW_POWERED_BY = True
NUM_SAMPLES = 0
NUM_VARIANTS = 'XYZ million'

MONGO = {
    'host': 'mongo',
    'port': 27017,
    'name': 'bravo'
}
DOWNLOAD_ALL_FILEPATH = ''
URL_PREFIX = ''

GOOGLE_ANALYTICS_TRACKING_ID = ''
SECRET_KEY = ''
GOOGLE_AUTH = False # True if app is using Google Auth 2.0
GOOGLE_LOGIN_CLIENT_ID = ''
GOOGLE_LOGIN_CLIENT_SECRET = ''
TERMS = True # True if app requires 'Terms of Use'. Can be used only if GOOGLE_AUTH is enabled.

EMAIL_WHITELIST = False # True if app has whitelisted emails. Can be used only if GOOGLE_AUTH is enabled.

API_GOOGLE_AUTH = False
API_IP_WHITELIST = [ '127.0.0.1' ]
API_VERSION = ''
API_DATASET_NAME = ''
API_COLLECTION_NAME = 'variants'
API_URL_PREFIX = '/api/' + API_VERSION
API_PAGE_SIZE = 1000
API_MAX_REGION = 250000
API_REQUESTS_RATE_LIMIT = ['1800/15 minute']

BRAVO_AUTH_SECRET = ''
BRAVO_ACCESS_SECRET = ''
BRAVO_AUTH_URL_PREFIX = '/api/' + API_VERSION + '/auth'

IGV_REFERENCE_PATH = '/data/genomes/'
IGV_CRAM_DIRECTORY = '/data/cram/'
IGV_CACHE_COLLECTION = 'igv_cache'
IGV_CACHE_DIRECTORY = '/data/cache/igv_cache/'
IGV_CACHE_LIMIT = 1000
BASE_COVERAGE_DIRECTORY = '/data/coverage/'
FASTA_URL = 'https://your.bravo.domain.name/genomes/hs38DH.fa'

ADMINS = [
    'email@email.email'
]

ADMIN = False # True if app is running in admin mode.
ADMIN_ALLOWED_IP = []

