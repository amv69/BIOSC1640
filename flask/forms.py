from flask_wtf import FlaskForm
from wtforms import IntegerField, SubmitField, StringField
from wtforms.validators import DataRequired, Optional

class Parameters(FlaskForm):
    fastA = StringField('FastA', validators=[DataRequired()])
    lengthOne = IntegerField('Length One', validators=[DataRequired()])
    lengthTwo = IntegerField('Length Two', validators=[DataRequired()])
    tempOne = IntegerField('Temperature One', validators=[DataRequired()])
    tempTwo = IntegerField('Temperature Two', validators=[DataRequired()])
    submit = SubmitField('SUBMIT')

