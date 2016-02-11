from flask import url_for, redirect, request
from rauth import OAuth2Service

import json, urllib2

class GoogleSignIn(object):
    def __init__(self, current_app):
        googleinfo = urllib2.urlopen('https://accounts.google.com/.well-known/openid-configuration')
        google_params = json.load(googleinfo)
        self.service = OAuth2Service(
            name='google',
            client_id=current_app.config['GOOGLE_LOGIN_CLIENT_ID'],
            client_secret=current_app.config['GOOGLE_LOGIN_CLIENT_SECRET'],
            authorize_url=google_params.get('authorization_endpoint'),
            base_url=google_params.get('userinfo_endpoint'),
            access_token_url=google_params.get('token_endpoint')
        )

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
            return None, None, None
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
