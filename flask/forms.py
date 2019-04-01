from flask_wtf import FlaskForm
from wtforms import IntegerField, SubmitField, StringField
from wtforms.validators import DataRequired, Optional
from wtforms.widgets import TextArea

class Parameters(FlaskForm):
    fastA = StringField('FastA', validators=[DataRequired()], widget=TextArea())
    lengthOne = StringField('Length Minimum ', validators=[DataRequired()])
    lengthTwo = StringField('Length Maximum ', validators=[DataRequired()])
    tempOne = StringField('Temperature Minimum ', validators=[DataRequired()])
    tempTwo = StringField('Temperature Maximum ', validators=[DataRequired()])
    maxTile = StringField('Maximum Tile Length ', validators=[DataRequired()])
    newSeq = StringField('Chosen Oligo')
    areaOne = StringField('Starting Position to Select')
    areaTwo = StringField('Ending Position')
    submit = SubmitField('SUBMIT')
    output = StringField('Output', widget=TextArea())

