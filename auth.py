import json
import urllib2

import requests
from flask import redirect, request, url_for
from rauth import OAuth2Service


class GoogleSignIn(object):
    def __init__(self, current_app):
        google_params = self._get_google_info()
        self.service = OAuth2Service(
            name = 'google',
            client_id = current_app.config['GOOGLE_LOGIN_CLIENT_ID'],
            client_secret = current_app.config['GOOGLE_LOGIN_CLIENT_SECRET'],
            authorize_url = google_params.get('authorization_endpoint'),
            base_url = google_params.get('userinfo_endpoint'),
            access_token_url = google_params.get('token_endpoint'))

    def _get_google_info(self):
        response = requests.get('https://accounts.google.com/.well-known/openid-configuration')
        response.raise_for_status()
        return json.loads(response.text)

    def authorize(self):
        return redirect(self.service.get_authorize_url(scope = 'email', response_type = 'code', prompt = 'select_account', redirect_uri = self.get_callback_url()))

    def get_callback_url(self):
        return url_for('.oauth_callback_google', _external = True)

    def callback(self):
        if 'code' not in request.args:
            return (None, None, None)
        data = {
            'code': request.args['code'],
            'grant_type': 'authorization_code',
            'redirect_uri': self.get_callback_url()
        }
        try:
            oauth_session = self.service.get_auth_session(data = data, decoder = json.loads)
        except:
            return (None, None, None)
        user = oauth_session.get('').json()
        return (user['name'], user['email'], user.get('picture'))
