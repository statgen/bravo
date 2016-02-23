import flask_wtf
from wtforms import TextAreaField, SubmitField, validators

class RequestStudyContactForm(flask_wtf.Form):
    info_requested = TextAreaField("Information you're interested in about people with this variant", [
        validators.Required("Please enter the information that you would like to receive about the people who have this variant.")
    ])
    submit = SubmitField("Send Request")
