from flask import url_for, redirect, request
from rauth import OAuth2Service

import json, urllib2
import requests

class GoogleSignIn(object):
    def __init__(self, current_app):
        google_params = self._get_google_info()
        self.service = OAuth2Service(
            name='google',
            client_id=current_app.config['GOOGLE_LOGIN_CLIENT_ID'],
            client_secret=current_app.config['GOOGLE_LOGIN_CLIENT_SECRET'],
            authorize_url=google_params.get('authorization_endpoint'),
            base_url=google_params.get('userinfo_endpoint'),
            access_token_url=google_params.get('token_endpoint')
        )

    def _get_google_info(self):
        # Previously I used:
        #  import urllib3.contrib.pyopenssl
        #  urllib3.contrib.pyopenssl.inject_into_urllib3()
        #  googleinfo = urllib2.urlopen('https://accounts.google.com/.well-known/openid-configuration')
        #  return json.load(googleinfo)
        r = requests.get('https://accounts.google.com/.well-known/openid-configuration')
        r.raise_for_status()
        return json.loads(r.text)

    def authorize(self):
        return redirect(self.service.get_authorize_url(
            scope='email',
            response_type='code',
            prompt='select_account',
            redirect_uri=self.get_callback_url())
            )

    def get_callback_url(self):
        return url_for('oauth_callback_google',
                        _external=True)

    def callback(self):
        if 'code' not in request.args:
            return (None, None)
        # The following two commands are based on `requests` and pass **kwargs to it.
        oauth_session = self.service.get_auth_session(
                data={'code': request.args['code'],
                      'grant_type': 'authorization_code',
                      'redirect_uri': self.get_callback_url()
                     },
                decoder = json.loads,
                verify = False
        )
        me = oauth_session.get('', verify=False).json()
        return (me['name'],
                me['email'])
