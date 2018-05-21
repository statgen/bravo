DEBUG = False
TESTING = False
PROXY = True # True if app is proxied by Apache or similar. 

BROWSER_NAME = 'Bravo'
DATASET_NAME = 'Example Dataset'
SHOW_POWERED_BY = True
NUM_SAMPLES = 0
NUM_VARIANTS = 'XYZ million'

MONGO = {
    'host': 'localhost',
    'port': 27017,
    'name': 'example'
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

BRAVO_API_VERSION = 'v1'
BRAVO_AUTH_SECRET = ''
BRAVO_ACCESS_SECRET = ''
BRAVO_AUTH_URL_PREFIX = '/api/' + BRAVO_API_VERSION + '/auth'
BRAVO_API_URL_PREFIX = '/api/' + BRAVO_API_VERSION

IGV_REFERENCE_PATH = ''
IGV_CRAM_DIRECTORY = ''
IGV_CACHE_COLLECTION = 'igv_cache'
IGV_CACHE_DIRECTORY = ''
IGV_CACHE_LIMIT = 1000

ADMINS = [
    'email@email.email'
]

ADMIN = False # True if app is running in admin mode.
ADMIN_ALLOWED_IP = []

#cov_dir = '/var/browser_coverage/topmed_freeze2_random1000_hg38/'
BASE_COVERAGE = []
#BASE_COVERAGE.extend({'bp-min-length':0,                  'path':path} for path in glob.glob(cov_dir + 'full/*.json.gz'))
#BASE_COVERAGE.extend({'bp-min-length':300, 'binned':True, 'path':path} for path in glob.glob(cov_dir + 'bin_25e-2/*.json.gz'))
#BASE_COVERAGE.extend({'bp-min-length':1000,'binned':True, 'path':path} for path in glob.glob(cov_dir + 'bin_50e-2/*.json.gz'))
#BASE_COVERAGE.extend({'bp-min-length':3000,'binned':True, 'path':path} for path in glob.glob(cov_dir + 'bin_75e-2/*.json.gz'))

