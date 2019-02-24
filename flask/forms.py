from flask_wtf import FlaskForm
from wtforms import IntegerField, SubmitField, StringField
from wtforms.validators import DataRequired, Optional

class Parameters(FlaskForm):
    fastA = StringField('fastA', validators=[Optional()])
    lengthOne = IntegerField('lengthOne', validators=[DataRequired()])
    lengthTwo = IntegerField('lengthTwo', validators=[DataRequired()])
    tempOne = IntegerField('tempOne', validators=[DataRequired()])
    tempTwo = IntegerField('tempTwo', validators=[DataRequired()])
    submit = SubmitField('submit')

