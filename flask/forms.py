from flask_wtf import FlaskForm
from wtforms import IntegerField, SubmitField, StringField
from wtforms.validators import DataRequired, Optional
from wtforms.widgets import TextArea

class Parameters(FlaskForm):
    fastA = StringField('FastA', validators=[DataRequired()], widget=TextArea())
    lengthOne = IntegerField('Length Minimum ', validators=[DataRequired()])
    lengthTwo = IntegerField('Length Maximum ', validators=[DataRequired()])
    tempOne = IntegerField('Temperature Minimum ', validators=[DataRequired()])
    tempTwo = IntegerField('Temperature Maximum ', validators=[DataRequired()])
    submit = SubmitField('SUBMIT')
    output = StringField('Output', widget=TextArea())

