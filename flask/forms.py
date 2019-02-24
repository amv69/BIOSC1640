from flask_wtf import FlaskForm
from wtforms import IntegerField, SubmitField, StringField
from wtforms.validators import DataRequired, Optional

class Parameters(FlaskForm):
    fastA = StringField('FastA', validators=[DataRequired()])
    lengthOne = IntegerField('Length Minimum', validators=[DataRequired()])
    lengthTwo = IntegerField('Length Maximum', validators=[DataRequired()])
    tempOne = IntegerField('Temperature Minimum', validators=[DataRequired()])
    tempTwo = IntegerField('Temperature Maximum', validators=[DataRequired()])
    submit = SubmitField('SUBMIT')

