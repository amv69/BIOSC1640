from flask import Flask, render_template
from forms import Parameters
import Oligo
app = Flask(__name__)

app.config['SECRET_KEY'] = '051259712ab12de981acb'

data = [
    {
        'fastA': 'ATCG',
        'lengthOne': '5',
        'lengthTwo': '10',
        'tempOne': '20',
        'tempTwo': '40',
    }

]

@app.route("/")
@app.route("/home")
def hello():
    return render_template('home.html', data=data)

@app.route("/about")
def about():
    return render_template('about.html', title='About')

@app.route("/submit", methods=['GET', 'POST'])
def submit():
    form = Parameters()
    if form.validate_on_submit()
        flash(f'Job Started!')
    return render_template('submit.html', title='Submit', form=form)

if __name__ == '__main__':
    app.run(debug=True)


